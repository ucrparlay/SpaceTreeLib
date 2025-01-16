#ifndef PSTP_BASE_TREE_IMPL_LEAF_OP_HPP_
#define PSTP_BASE_TREE_IMPL_LEAF_OP_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace pstp {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::ExtractPointsInLeaf(
    Node* T, Range Out) {
  Leaf* TL = static_cast<Leaf*>(T);
  if (TL->is_dummy) {
    std::ranges::fill_n(Out.begin(), TL->size, TL->pts[0]);
  } else {
    std::ranges::copy(TL->pts.begin(), TL->pts.begin() + TL->size, Out.begin());
  }
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::InsertPoints2Leaf(
    Node* T, Slice In) {
  Leaf* TL = static_cast<Leaf*>(T);
  if (TL->is_dummy) {
    T->size += In.size();
    return T;
  }

  if (TL->pts.size() == 0) {
    TL->pts = Points::uninitialized(kLeaveWrap);
  }
  for (size_t i = 0; i < In.size(); i++) {
    TL->pts[TL->size + i] = In[i];
  }
  TL->size += In.size();
  return T;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename RT>
RT BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DeletePoints4Leaf(
    Node* T, Slice In) {
  assert(T->size >= In.size());
  Leaf* TL = static_cast<Leaf*>(T);

  if (TL->is_dummy) {
    assert(In.size() <= T->size);
    TL->size -= In.size();  // WARN: this assumes that In\in T
    if (TL->size == 0) {
      TL->is_dummy = false;
      TL->pts = Points::uninitialized(kLeaveWrap);
    }

    if constexpr (std::same_as<RT, Node*>) {
      return T;
    } else if constexpr (std::same_as<RT, NodeBox>) {
      return NodeBox(T, T->size ? Box(TL->pts[0], TL->pts[0]) : GetEmptyBox());
    } else {
      static_assert(std::same_as<RT, Node*> || std::same_as<RT, NodeBox>);
    }
  }

  auto it = TL->pts.begin(), end = TL->pts.begin() + TL->size;
  for (size_t i = 0; i < In.size(); i++) {
    it = std::ranges::find(TL->pts.begin(), end, In[i]);
    assert(it != end);
    std::ranges::iter_swap(it, --end);
  }

  assert(std::cmp_equal(std::distance(TL->pts.begin(), end),
                        TL->size - In.size()));
  TL->size -= In.size();
  assert(TL->size >= 0);

  if constexpr (std::same_as<RT, Node*>) {
    return T;
  } else if constexpr (std::same_as<RT, NodeBox>) {
    return NodeBox(T, GetBox(TL->pts.cut(0, TL->size)));
  } else {
    static_assert(std::same_as<RT, Node*> || std::same_as<RT, NodeBox>);
  }
}

// NOTE: diff points from the leaf using std::set_difference
// {1, 2, 5, 5, 5, 9} ∖ {2, 5, 7} == {1, 5, 5, 9}
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename RT>
RT BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DiffPoints4Leaf(
    Node* T, Slice In) {
  Leaf* TL = static_cast<Leaf*>(T);

  if (TL->is_dummy) {
    size_t cnt = parlay::count(In, TL->pts[0]);
    TL->size = TL->size >= cnt ? TL->size - cnt : 0;
    assert(TL->size >= 0);
    if (TL->size == 0) {  // set points to normal leaf when all points in dummy
                          // node has been deleted
      TL->is_dummy = false;
      TL->pts = Points::uninitialized(kLeaveWrap);
    }

    if constexpr (std::same_as<RT, Node*>) {
      return T;
    } else if constexpr (std::same_as<RT, NodeBox>) {
      return NodeBox(T, TL->size ? Box(TL->pts[0], TL->pts[0]) : GetEmptyBox());
    } else {
      static_assert(std::same_as<RT, Node*> || std::same_as<RT, NodeBox>);
    }
  }

  // NOTE: for normal leaf, need to check whether all Points are in the leaf
  auto diff_res = std::ranges::set_difference(
      parlay::sort(TL->pts.cut(0, TL->size)), parlay::sort(In), TL->pts.begin(),
      [](Point const& p1, Point const& p2) { return p1 < p2; });
  TL->size = std::ranges::distance(TL->pts.begin(), diff_res.out);

  if constexpr (std::same_as<RT, Node*>) {
    return T;
  } else if constexpr (std::same_as<RT, NodeBox>) {
    return NodeBox(T, GetBox(TL->pts.cut(0, TL->size)));
  } else {
    static_assert(std::same_as<RT, Node*> || std::same_as<RT, NodeBox>);
  }
}
}  // namespace pstp

#endif  // PSTP_BASE_TREE_IMPL_LEAF_OP_HPP_
