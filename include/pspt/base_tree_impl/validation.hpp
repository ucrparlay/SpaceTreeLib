#ifndef PSPT_BASE_TREE_IMPL_VALIDATION_HPP_
#define PSPT_BASE_TREE_IMPL_VALIDATION_HPP_

#include <parlay/parallel.h>

#include <concepts>
#include <cwchar>

#include "../base_tree.h"

namespace pspt {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckBox(Node* T,
                                                              Box const& box) {
  if (T->is_leaf) {
    [[maybe_unused]] auto const* TL = static_cast<Leaf const*>(T);
    [[maybe_unused]] auto const node_box = GetBox<Leaf, Interior>(T);
    if constexpr (std::same_as<decltype(TL->GetSplit()),
                               Box>) {  // whether storing a bounding box
      SameBox(node_box, box);
    } else {
      WithinBox(node_box, box);
    }
    return node_box;
  }
  Interior* TI = static_cast<Interior*>(T);
  assert(!TI->GetParallelFlagIniStatus());  // NOTE: ensure that uninitialized
                                            // force parallelism
  if constexpr (IsBinaryNode<Interior> &&
                std::same_as<typename Interior::ST,
                             HyperPlane>) {  // NOTE: use hyperplane as splitter
    Box lbox(box), rbox(box);
    lbox.second.pnt[TI->split.second] = TI->split.first;
    rbox.first.pnt[TI->split.second] = TI->split.first;
    Box const left_return_box = CheckBox<Leaf, Interior>(TI->left, lbox);
    Box const right_return_box = CheckBox<Leaf, Interior>(TI->right, rbox);
    Box const new_box = GetBox(left_return_box, right_return_box);

    assert(WithinBox(left_return_box, lbox));
    assert(WithinBox(right_return_box, rbox));
    assert(WithinBox(new_box, box));
    return new_box;
  } else if constexpr (IsBinaryNode<Interior> &&
                       std::same_as<typename Interior::ST,
                                    Box>) {  // use box as splitter
    Box const left_return_box =
        CheckBox<Leaf, Interior>(TI->left, GetSplit<Leaf, Interior>(TI->left));
    Box const right_return_box = CheckBox<Leaf, Interior>(
        TI->right, GetSplit<Leaf, Interior>(TI->right));
    Box const new_box = GetBox(left_return_box, right_return_box);
    assert(SameBox(left_return_box, GetSplit<Leaf, Interior>(TI->left)));
    assert(SameBox(right_return_box, GetSplit<Leaf, Interior>(TI->right)));
    assert(SameBox(new_box, box));
    return new_box;
  } else if (IsMultiNode<Interior>) {
    BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));
    BoxSeq return_box_seq(new_box.size());
    assert(new_box.size() == TI->tree_nodes.size());
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      return_box_seq[i] =
          CheckBox<Leaf, Interior>(TI->tree_nodes[i], new_box[i]);
      assert(WithinBox(return_box_seq[i], new_box[i]));
    }
    auto return_box = GetBox(return_box_seq);
    assert(WithinBox(return_box, box));
    return return_box;
  }
}

// NOTE: 1: whether the subtree is within the circle
// 2: separation of nodes per level
// 3: decrease level per level
// 4: nesting
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename Interior::CircleType
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckCover(
    Node* T, typename Interior::CircleType const& circle) {
  using CircleType = typename Interior::CircleType;

  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    auto leaf_pnt_circle = GetCircle<CircleType>(TL->pts.cut(0, TL->size));
    assert(CircleWithinCircle(leaf_pnt_circle, circle));
    return leaf_pnt_circle;
  }

  auto TI = static_cast<Interior*>(T);
  assert(circle == TI->GetCoverCircle());
  parlay::sequence<CircleType> return_circle_seq(TI->split.size());
  for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
    auto next_circle = CircleType{
        TI->split[i], static_cast<DepthType>(
                          static_cast<DepthType>(TI->GetCoverCircle().level) -
                          static_cast<DepthType>(1))};
    return_circle_seq[i] =
        CheckCover<Leaf, Interior>(TI->tree_nodes[i], next_circle);
    assert(CircleWithinCircle(return_circle_seq[i], next_circle));
  }
  auto const return_circle = GetCircle<CircleType>(return_circle_seq);
  // static_assert(std::same_as<decltype(return_circle), decltype(circle)>);
  assert(CircleWithinCircle(return_circle, circle));
  return return_circle;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckSize(Node* T) {
  if (T->is_leaf) {
    assert(static_cast<Leaf*>(T)->is_dummy ||
           static_cast<Leaf*>(T)->size <= kLeaveWrap);
    return T->size;
  }
  if constexpr (IsBinaryNode<Interior>) {
    Interior* TI = static_cast<Interior*>(T);
    assert(!TI->GetParallelFlagIniStatus());
    size_t l, r;
    parlay::par_do([&l, &TI] { l = CheckSize<Leaf, Interior>(TI->left); },
                   [&r, &TI] { r = CheckSize<Leaf, Interior>(TI->right); });
    assert(l + r == T->size);
    return T->size;
  } else {
    // assert(IsMultiNode<Interior>);
    Interior* TI = static_cast<Interior*>(T);
    assert(!TI->GetParallelFlagIniStatus());
    size_t sum = 0;
    for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
      sum += CheckSize<Leaf, Interior>(TI->tree_nodes[i]);
    }
    assert(sum == T->size);
    return T->size;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
void BaseTree<Point, DerivedTree, kSkHeight,
              kImbaRatio>::CheckTreeSameSequential(Node* T, int dim) {
  if (T->is_leaf) {
    // assert( PickRebuildDim( T, kDim ) == dim );
    return;
  }
  if constexpr (IsBinaryNode<Interior>) {
    Interior* TI = static_cast<Interior*>(T);
    if (TI->split.second != dim) {
      std::cout << int(TI->split.second) << " " << int(dim) << " " << TI->size;
    }
    assert(TI->split.second == dim);
    // TODO: maybe need to add the split rule in the base tree?
    dim = (dim + 1) % kDim;
    parlay::par_do_if(
        T->size > 1000,
        [&]() { CheckTreeSameSequential<Leaf, Interior>(TI->left, dim); },
        [&]() { CheckTreeSameSequential<Leaf, Interior>(TI->right, dim); });
  } else {
    assert(IsMultiNode<Interior>);
    Interior* TI = static_cast<Interior*>(T);
    assert(std::cmp_equal((1 << TI->split.size()), TI->tree_nodes.size()));
    for (size_t i = 0; i < TI->split.size(); i++) {
      assert(TI->split[i].second == dim);
      dim += 1;
    }
    assert(dim == kDim);
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      CheckTreeSameSequential<Leaf, Interior>(TI->tree_nodes[i], 0);
    }
  }
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, typename SplitRule>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Validate() {
  std::cout << ">>> begin validate tree\n" << std::flush;
  if constexpr (IsBinaryNode<Interior> || IsMultiNode<Interior>) {
    if (LegalBox(CheckBox<Leaf, Interior>(this->root_, this->tree_box_))) {
      std::cout << "Correct bounding Box\n" << std::flush;
    } else {
      std::cout << "wrong bounding Box\n" << std::flush;
      abort();
    }

    // NOTE: used to check rotate dimension
    // For kdtree binary node, the dummy node may break the rotation manner,
    // since if current dimension is un-splitable, one has to switch to another
    // dimension
    if constexpr (IsRotateDimSplit<typename SplitRule::DimRuleType> &&
                  IsMultiNode<Interior>) {
      CheckTreeSameSequential<Leaf, Interior>(this->root_, 0);
      std::cout << "Correct rotate dimension\n" << std::flush;
    }
  } else if constexpr (IsDynamicNode<Interior>) {
    auto root_cover_circle =
        this->root_->is_leaf
            ? GetCircle<typename Interior::CircleType>(
                  static_cast<Leaf*>(this->root_)
                      ->pts.cut(0, this->root_->size))
            : static_cast<Interior*>(this->root_)->GetCoverCircle();

    static_assert(std::same_as<decltype(root_cover_circle),
                               typename Interior::CircleType>);

    if (LegalCircle(
            CheckCover<Leaf, Interior>(this->root_, root_cover_circle))) {
      std::cout << "Correct cover circle\n" << std::flush;
    } else {
      std::cout << "wrong bounding Box\n" << std::flush;
      abort();
    }
  } else {
    assert(false);
  }

  if (CheckSize<Leaf, Interior>(this->root_) == this->root_->size) {
    std::cout << "Correct size\n" << std::flush;
  } else {
    std::cout << "wrong tree size\n" << std::flush;
    abort();
  }
  std::cout << "<<< end validate tree\n" << std::flush;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetTreeHeight() {
  size_t deep = 0;
  return GetMaxTreeDepth<Leaf, Interior>(this->root_, deep);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetMaxTreeDepth(
    Node* T, size_t deep) {
  if (T->is_leaf) {
    return deep;
  }

  Interior* TI = static_cast<Interior*>(T);
  if constexpr (IsBinaryNode<Interior>) {
    int l = GetMaxTreeDepth<Leaf, Interior>(TI->left, deep + 1);
    int r = GetMaxTreeDepth<Leaf, Interior>(TI->right, deep + 1);
    return std::max(l, r);
  } else if constexpr (IsMultiNode<Interior>) {
    size_t max_depth = 0;
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      max_depth = std::max(max_depth, GetMaxTreeDepth<Leaf, Interior>(
                                          TI->tree_nodes[i], deep + 1));
    }
    return max_depth;
  } else {
    return 0;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
double BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetAveTreeHeight() {
  // std::cout << IsBinaryNode<Interior> << std::endl;
  parlay::sequence<size_t> heights(this->root_->size);
  size_t idx = 0;
  CountTreeHeights<Leaf, Interior>(this->root_, 0, idx, heights);
  // auto kv = parlay::histogram_by_key( heights.cut( 0, idx ) );
  // std::sort( kv.begin(), kv.end(),
  //            [&]( auto a, auto b ) { return a.first < b.first; } );
  // for ( auto i : kv )
  //     std::cout << i.first << " " << i.second << std::endl;
  return double(1.0 * parlay::reduce(heights.cut(0, idx)) / idx);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CountTreeNodesNum(
    Node* T) {
  if (T->is_leaf) {
    return 1;
  }

  Interior* TI = static_cast<Interior*>(T);
  if constexpr (IsBinaryNode<Interior>) {
    size_t l, r;
    parlay::par_do([&]() { l = CountTreeNodesNum<Leaf, Interior>(TI->left); },
                   [&]() { r = CountTreeNodesNum<Leaf, Interior>(TI->right); });
    return l + r + 1;
  } else {
    size_t sum = 0;
    for (int i = 0; i < TI->tree_nodes.size(); i++) {
      sum += CountTreeNodesNum<Leaf, Interior>(TI->tree_nodes[i]);
    }
    return sum + 1;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CountTreeHeights(
    Node* T, size_t deep, size_t& idx, parlay::sequence<size_t>& heights) {
  if (T->is_leaf) {
    heights[idx++] = deep;
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  if constexpr (IsBinaryNode<Interior>) {
    CountTreeHeights<Leaf, Interior>(TI->left, deep + 1, idx, heights);
    CountTreeHeights<Leaf, Interior>(TI->right, deep + 1, idx, heights);
  } else if constexpr (IsMultiNode<Interior>) {
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      CountTreeHeights<Leaf, Interior>(TI->tree_nodes[i], deep + 1, idx,
                                       heights);
    }
  } else {
    return;
  }
  return;
}

}  // namespace pspt

#endif  // PSPT_BASE_TREE_IMPL_VALIDATION_HPP_
