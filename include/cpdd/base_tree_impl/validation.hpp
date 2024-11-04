#pragma once

#include <parlay/parallel.h>

#include <cwchar>

#include "../base_tree.h"

namespace cpdd {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::CheckBox(Node* T,
                                                              Box const& box) {
  if (T->is_leaf) {
    return GetBox<Leaf, Interior>(T);
  }
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (IsBinaryNode<Interior>) {
    Box lbox(box), rbox(box);
    lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
    rbox.first.pnt[TI->split.second] = TI->split.first;
    Box const left_return_box = CheckBox<Leaf, Interior>(TI->left, lbox);
    Box const right_return_box = CheckBox<Leaf, Interior>(TI->right, rbox);
    Box const new_box = GetBox(left_return_box, right_return_box);

    assert(WithinBox(left_return_box, lbox));
    assert(WithinBox(right_return_box, rbox));
    assert(WithinBox(new_box, box));
    return new_box;
  } else {
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
    assert(IsMultiNode<Interior>);
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
      std::cout << int(TI->split.second) << " " << int(dim) << " " << TI->size
                << std::endl;
    }
    assert(TI->split.second == dim);
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
  std::cout << ">>> begin validate tree" << std::endl << std::flush;
  if (LegalBox(CheckBox<Leaf, Interior>(this->root_, this->tree_box_))) {
    std::cout << "Correct bounding Box" << std::endl << std::flush;
  } else {
    std::cout << "wrong bounding Box" << std::endl << std::flush;
    abort();
  }

  // NOTE: used to check rotate dimension
  // For kdtree binary node, the dummy node may break the rotation manner,
  // since if current dimension is un-splitable, one has to switch to another
  // dimension
  if constexpr (IsRotateDimSplit<SplitRule> && IsMultiNode<Interior>) {
    CheckTreeSameSequential<Leaf, Interior>(this->root_, 0);
    std::cout << "Correct rotate dimension" << std::endl << std::flush;
  }

  if (CheckSize<Leaf, Interior>(this->root_) == this->root_->size) {
    std::cout << "Correct size" << std::endl << std::flush;
  } else {
    std::cout << "wrong tree size" << std::endl << std::flush;
    abort();
  }
  std::cout << "<<< end validate tree" << std::endl << std::flush;
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
  } else {
    size_t max_depth = 0;
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      max_depth = std::max(max_depth, GetMaxTreeDepth<Leaf, Interior>(
                                          TI->tree_nodes[i], deep + 1));
    }
    return max_depth;
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
  } else {
    for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
      CountTreeHeights<Leaf, Interior>(TI->tree_nodes[i], deep + 1, idx,
                                       heights);
    }
  }
  return;
}

}  // namespace cpdd
