#ifndef PSI_BASE_TREE_IMPL_VALIDATION_HPP_
#define PSI_BASE_TREE_IMPL_VALIDATION_HPP_

#include <parlay/parallel.h>

#include <concepts>
#include <cwchar>

#include "../base_tree.h"
#include "parlay/primitives.h"

namespace psi {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckBox(Node* T,
                                                              Box const& box) {
  if (T->is_leaf) {
    auto const* TL = static_cast<Leaf const*>(T);
    auto const node_box = GetBox<Leaf, Interior>(T);
    if constexpr (HasBox<typename Leaf::AT>) {  // whether has bb
      assert(SameBox(node_box, TL->GetBox()));
    } else {
      assert(WithinBox(node_box, box));
    }
    return node_box;
  }
  Interior* TI = static_cast<Interior*>(T);
  assert(!TI->GetParallelFlagIniStatus());  // NOTE: ensure that uninitialized
                                            // force parallelism
  if constexpr (IsBinaryNode<Interior> &&
                !HasBox<typename Interior::AT>) {  // NOTE: use hyperplane as
                                                   // splitter
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
                       HasBox<typename Interior::AT>) {  // kd with box
    // std::cout << " has box " << std::endl;
    auto left_box = RetriveBox<Leaf, Interior>(TI->left);
    auto right_box = RetriveBox<Leaf, Interior>(TI->right);
    Box const left_return_box = CheckBox<Leaf, Interior>(TI->left, left_box);
    Box const right_return_box = CheckBox<Leaf, Interior>(TI->right, right_box);
    Box const new_box = GetBox(left_return_box, right_return_box);
    assert(SameBox(left_return_box, RetriveBox<Leaf, Interior>(TI->left)));
    assert(SameBox(right_return_box, RetriveBox<Leaf, Interior>(TI->right)));
    assert(SameBox(new_box, TI->GetBox()));
    return new_box;
  } else if (IsMultiNode<Interior> &&
             !HasBox<typename Interior::AT>) {  // orth without box
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
  } else if (IsMultiNode<Interior> &&
             HasBox<typename Interior::AT>) {  // orth with box
    BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));
    BoxSeq return_box_seq(new_box.size());
    assert(new_box.size() == TI->tree_nodes.size());
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      return_box_seq[i] =
          CheckBox<Leaf, Interior>(TI->tree_nodes[i], new_box[i]);
      assert(SameBox(return_box_seq[i],
                     RetriveBox<Leaf, Interior>(TI->tree_nodes[i])));
    }
    auto return_box = GetBox(return_box_seq);
    assert(SameBox(return_box, TI->GetBox()));
    return return_box;
  } else {
    assert(false);
  }
}

// NOTE: 1: whether the subtree is within the circle
// 2: separation of nodes per level
// 3: decrease level per level
// 4: nesting
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Points
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckCover(
    Node* T, typename Interior::CircleType const& level_cover_circle) {
  using TreeCircle = typename Interior::CircleType;

  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    for (auto const& p : TL->GetPoints()) {
      assert(WithinCircle(p, level_cover_circle));
    }
    return TL->GetPoints();
  }

  auto TI = static_cast<Interior*>(T);
  assert(level_cover_circle == TI->GetCoverCircle());

  // NOTE: nesting
  assert(std::ranges::count(TI->split, level_cover_circle.center) == 1);

  // NOTE: separation
  assert(std::all_of(
      TI->split.begin(), TI->split.end(), [&](auto const& center_i) {
        auto i = &center_i - &TI->split[0];
        auto tmp_circle = CoverCircle{center_i, level_cover_circle.level - 1};
        return std::all_of(TI->split.begin() + i + 1, TI->split.end(),
                           [&](auto const& center_j) {
                             // PERF: should be Gt as a point in a boundary
                             // should falls within that circle
                             return Num::Gt(
                                 P2PDistanceSquare(center_i, center_j),
                                 tmp_circle.GetRadiusSquare());
                           });
      }));

  // NOTE: covering
  // the cover circle of all points within the subtree should within the cover
  // circle of the node
  parlay::sequence<Points> return_points_seq(TI->tree_nodes.size());
  for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
    auto next_circle = TreeCircle{
        TI->split[i], static_cast<DepthType>(TI->GetCoverCircle().level - 1)};
    return_points_seq[i] =
        CheckCover<Leaf, Interior>(TI->tree_nodes[i], next_circle);
    parlay::all_of(return_points_seq[i], [&](auto const& p) {
      return WithinCircle(p, level_cover_circle);
    });
    // for (auto const& p : return_points_seq[i]) {
    //   assert(WithinCircle(p, next_circle));
    // }
  }
  return parlay::flatten(return_points_seq);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckSize(Node* T) {
  if (T->is_leaf) {
    // assert(static_cast<Leaf*>(T)->is_dummy ||
    //        static_cast<Leaf*>(T)->size <= kLeaveWrap);
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

  // check size
  if (CheckSize<Leaf, Interior>(this->root_) == this->root_->size) {
    std::cout << "Correct size\n" << std::flush;
  } else {
    std::cout << "wrong tree size\n" << std::flush;
    abort();
  }

  // tree property
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
    // BUG: below is wrong, the leaf cover circle should centered at one node.
    // Cannot compute it, but to retrive from the tree
    auto root_cover_circle =
        this->root_->is_leaf
            ? GetCircle<typename Interior::CircleType>(
                  static_cast<Leaf*>(this->root_)
                      ->pts.cut(0, this->root_->size))
            : static_cast<Interior*>(this->root_)->GetCoverCircle();

    static_assert(std::same_as<decltype(root_cover_circle),
                               typename Interior::CircleType>);

    if (parlay::all_of(
            CheckCover<Leaf, Interior>(this->root_, root_cover_circle),
            [&](auto const& p) {
              return WithinCircle(p, root_cover_circle);
            })) {
      std::cout << "Correct cover circle\n" << std::flush;
    } else {
      std::cout << "wrong bounding Box\n" << std::flush;
      abort();
    }
  } else {
    assert(false);
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

}  // namespace psi

#endif  // PSI_BASE_TREE_IMPL_VALIDATION_HPP_
