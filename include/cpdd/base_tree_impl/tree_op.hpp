#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../base_tree.h"
#include "cpdd/dependence/concepts.h"
#include "cpdd/dependence/tree_node.h"
#include "parlay/primitives.h"
#include "parlay/slice.h"

namespace cpdd {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildWithInsert(
    Node* T, Slice In, Args&&... args) {
  Points wx, wo;
  PrepareRebuild<Leaf, Interior>(T, In, wx, wo);
  static_assert(
      std::is_invocable_v<decltype(&DerivedTree::BuildRecursive), DerivedTree*,
                          Slice, Slice, Args&&..., Box>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(wx), parlay::make_slice(wo),
      std::forward<Args>(args)..., GetBox(parlay::make_slice(wx)));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, bool granularity, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildSingleTree(
    Node* T, Args&&... args) {
  Points wx, wo;
  PrepareRebuild<Leaf, Interior, granularity>(T, wx, wo);
  static_assert(std::is_invocable_v<decltype(&DerivedTree::BuildRecursive),
                                    DerivedTree*, Slice, Slice, Args&&...>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(wx), parlay::make_slice(wo),
      std::forward<Args>(args)...);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildTreeRecursive(
    Node* T, PrepareFunc&& prepare_func, Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  if (ImbalanceNode(TI->left->size, TI->size)) {
    return RebuildSingleTree<Leaf, Interior, granularity>(
        T, std::forward<Args>(args)...);
  }

  auto const [left_args, right_args] =
      prepare_func(T, std::forward<Args>(args)...);

  Node *L, *R;
  parlay::par_do_if(
      ForceParallelRecursion<Interior, granularity>(TI),
      [&] {
        L = std::apply(
            [&](auto&&... left_args) {
              return RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->left, prepare_func,
                  std::forward<decltype(left_args)>(left_args)...);
            },
            left_args);
      },
      [&] {
        R = std::apply(
            [&](auto&&... right_args) {
              return RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->right, prepare_func,
                  std::forward<decltype(right_args)>(right_args)...);
            },
            right_args);
      });

  TI->ResetParallelFlag();
  UpdateInterior<Interior>(T, L, R);
  return T;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildTreeRecursive(
    Node* T, PrepareFunc&& prepare_func, Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  if (SparcyNode(0,
                 TI->size)) {  // NOTE: after diff points from the tree, all
                               // points has been removed, so there the remove
                               // points are 0, simply use it to check whether
                               // the number of points are below kLeaveWrap
    return RebuildSingleTree<Leaf, Interior, granularity>(
        T, std::forward<Args>(args)...);
  }

  typename Interior::NodeArr new_nodes;
  if (ForceParallelRecursion<Interior, granularity>(TI)) {
    parlay::parallel_for(0, TI->tree_nodes.size(), [&](BucketType i) {
      auto const new_args = prepare_func(T, i, std::forward<Args>(args)...);
      std::apply(
          [&](auto&&... new_args) {
            new_nodes[i] = RebuildTreeRecursive<Leaf, Interior, granularity>(
                TI->tree_nodes[i], prepare_func,
                std::forward<decltype(new_args)>(new_args)...);
          },
          new_args);
    });
  } else {
    for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
      auto const new_args = prepare_func(T, i, std::forward<Args>(args)...);
      std::apply(
          [&](auto&&... new_args) {
            new_nodes[i] = RebuildTreeRecursive<Leaf, Interior, granularity>(
                TI->tree_nodes[i], prepare_func,
                std::forward<decltype(new_args)>(new_args)...);
          },
          new_args);
    }
  }

  TI->ResetParallelFlag();
  UpdateInterior<Interior>(T, new_nodes);
  return T;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <SupportsForceParallel Interior, bool granularity>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::ForceParallelRecursion(Interior const* TI) {
  return (granularity && TI->size > kSerialBuildCutoff) ||
         (!granularity && TI->ForceParallel());
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode Interior, typename Base>
Base BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<Base> const& tree_nodes) {
  if (idx > kPivotNum) {
    assert(idx - kPivotNum - 1 < kBucketNum);
    return tree_nodes[idx - kPivotNum - 1];
  }

  Base L, R;
  L = BuildInnerTree<Interior>(idx << 1, pivots, tree_nodes);
  R = BuildInnerTree<Interior>(idx << 1 | 1, pivots, tree_nodes);

  if constexpr (IsPointerToNode<Base>) {
    return AllocInteriorNode<Interior>(L, R, pivots[idx],
                                       typename Interior::AT());
  } else if constexpr (IsNodeBox<Base, Point>) {
    auto new_box = GetBox(L.second, R.second);
    return Base(AllocInteriorNode<Interior>(L.first, R.first, new_box,
                                            typename Interior::AT()),
                new_box);
  } else {
    static_assert(IsPointerToNode<Base> || IsNodeBox<Base, Point>);
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsMultiNode Interior>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<Node*> const& tree_nodes) {
  assert(idx < kPivotNum + kBucketNum + 1);

  if (idx > kPivotNum) {
    return tree_nodes[idx - kPivotNum - 1];
  }

  typename DerivedTree::OrthNodeArr multi_nodes;
  typename DerivedTree::Splitter split;
  for (DimsType i = 0; i < DerivedTree::kNodeRegions; ++i) {
    multi_nodes[i] = BuildInnerTree<Interior>(
        idx * DerivedTree::kNodeRegions + i, pivots, tree_nodes);
  }
  for (DimsType i = 0; i < DerivedTree::kSplitterNum; ++i) {
    split[i] = pivots[idx * (1 << i)];
    assert(i == 0 || pivots[idx * (1 << i)] == pivots[idx * (1 << i) + 1]);
  }

  return AllocInteriorNode<Interior>(multi_nodes, split,
                                     typename Interior::AT());
}

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
template <typename Leaf, IsBinaryNode Interior, typename Range,
          bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::FlattenRec(
    Node* T, Range Out) {
  assert(T->size == Out.size());

  if (T->size == 0) return;

  if (T->is_leaf) {
    ExtractPointsInLeaf<Leaf>(T, Out);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  assert(TI->size == TI->left->size + TI->right->size);
  parlay::par_do_if(
      // WARN: check parallelisim using node size can be biased
      ForceParallelRecursion<Interior, granularity>(TI),
      [&]() {
        FlattenRec<Leaf, Interior>(TI->left, Out.cut(0, TI->left->size));
      },
      [&]() {
        FlattenRec<Leaf, Interior>(TI->right,
                                   Out.cut(TI->left->size, TI->size));
      });

  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::FlattenRec(
    Node* T, Range Out) {
  if (T->size != Out.size()) {
    std::cout << "T->size: " << T->size << " Out.size(): " << Out.size()
              << std::endl;
  }
  assert(T->size == Out.size());

  if (T->size == 0) return;

  if (T->is_leaf) {
    ExtractPointsInLeaf<Leaf>(T, Out);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);

  assert(TI->size ==
         std::accumulate(
             TI->tree_nodes.begin(), TI->tree_nodes.end(),
             static_cast<size_t>(0),
             [](size_t acc, Node* n) -> size_t { return acc + n->size; }));

  if (ForceParallelRecursion<Interior, granularity>(TI)) {
    parlay::parallel_for(0, TI->tree_nodes.size(), [&](BucketType i) {
      size_t start = 0;
      for (BucketType j = 0; j < i; ++j) {
        start += TI->tree_nodes[j]->size;
      }
      FlattenRec<Leaf, Interior, Range>(
          TI->tree_nodes[i], Out.cut(start, start + TI->tree_nodes[i]->size));
    });
  } else {
    size_t start = 0;
    for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
      FlattenRec<Leaf, Interior, Range>(
          TI->tree_nodes[i], Out.cut(start, start + TI->tree_nodes[i]->size));
      start += TI->tree_nodes[i]->size;
    }
  }

  return;
}

// NOTE: for multi node @T, it only flatten the subtree with id @idx to @Out
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::PartialFlatten(
    Node* T, Range Out, BucketType idx) {
  if (idx == 1) {
    assert(T->size == Out.size());
    FlattenRec<Leaf, Interior>(T, Out.cut(0, T->size));
    return;
  } else if (idx >= Interior::kRegions) {
    Node* ns = static_cast<Interior*>(T)->tree_nodes[idx - Interior::kRegions];
    assert(ns->size == Out.size());
    FlattenRec<Leaf, Interior>(ns, Out.cut(0, ns->size));
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  size_t l_size = TI->MergeSize(idx << 1), r_size = TI->MergeSize(idx << 1 | 1);
  assert(l_size + r_size == Out.size());
  parlay::par_do_if(
      ForceParallelRecursion<Interior, granularity>(static_cast<Interior*>(T)),
      [&]() {
        PartialFlatten<Leaf, Interior>(T, Out.cut(0, l_size), idx << 1);
      },
      [&]() {
        PartialFlatten<Leaf, Interior>(T, Out.cut(l_size, l_size + r_size),
                                       idx << 1 | 1);
      });
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::ExpandMultiNode(
    typename MN::ST const& split, BucketType idx, BucketType deep,
    parlay::sequence<Node*> const& tree_nodes) {
  if (idx >= MN::kRegions) {
    return tree_nodes[idx - MN::kRegions];
  }
  Node *L, *R;
  L = ExpandMultiNode<BN, MN>(split, idx * 2, deep + 1, tree_nodes);
  R = ExpandMultiNode<BN, MN>(split, idx * 2 + 1, deep + 1, tree_nodes);
  auto o = AllocInteriorNode<BN>(L, R, split[deep], typename BN::AT());
  assert(o->size == L->size + R->size);
  return o;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode BN, IsMultiNode MN>
  requires std::same_as<typename BN::ST, typename MN::ST::value_type>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Expand2Binary(
    Node* T) {
  if (T->is_leaf) {
    return T;
  }
  MN* MI = static_cast<MN*>(T);
  parlay::sequence<Node*> tree_nodes(MN::kRegions);
  assert(MI->tree_nodes.size() == MN::kRegions);
  for (BucketType i = 0; i < MN::kRegions; ++i) {
    tree_nodes[i] = Expand2Binary<BN, MN>(MI->tree_nodes[i]);
  }
  assert(std::accumulate(tree_nodes.begin(), tree_nodes.end(), 0,
                         [](size_t acc, Node* n) -> size_t {
                           return acc + n->size;
                         }) == MI->size);
  auto split = MI->split;
  // FreeNode<MN>(T);
  return ExpandMultiNode<BN, MN>(split, 1, 0, tree_nodes);
}

// NOTE: update the info of T by new children L and R
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode Interior>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, Node* L, Node* R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  TI->ResetParallelFlag();
  TI->size = L->size + R->size;
  TI->left = L;
  TI->right = R;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode Interior>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, NodeBox const& L, NodeBox const& R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  TI->ResetParallelFlag();
  TI->size = L.first->size + R.first->size;
  TI->left = L.first;
  TI->right = R.first;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsMultiNode Interior>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, typename Interior::NodeArr const& new_nodes) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  TI->ResetParallelFlag();
  TI->size = std::accumulate(
      new_nodes.begin(), new_nodes.end(), 0,
      [](size_t acc, Node* n) -> size_t { return acc + n->size; });
  TI->tree_nodes = new_nodes;
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
    if (TL->size == 0) {
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

  // NOTE: need to check whether all Points are in the Leaf
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

}  // namespace cpdd
