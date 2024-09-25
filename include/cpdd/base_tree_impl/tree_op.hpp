
#pragma once

#include <cassert>
#include <numeric>
#include <utility>
#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"
#include "parlay/slice.h"

namespace cpdd {
template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, typename Interior, typename... Args>
Node* BaseTree<Point, DerivedTree, kBDO>::RebuildWithInsert(Node* T, Slice In,
                                                            Args&&... args) {
    Points wx, wo;
    PrepareRebuild<Leaf, Interior>(T, In, wx, wo);
    static_assert(
        std::is_invocable_v<decltype(&DerivedTree::BuildRecursive),
                            DerivedTree*, Slice, Slice, Args&&..., Box>);
    return static_cast<DerivedTree*>(this)->BuildRecursive(
        parlay::make_slice(wx), parlay::make_slice(wo),
        std::forward<Args>(args)..., GetBox(parlay::make_slice(wx)));
}

// TODO: if the bounding box has already been provided, we should not return one
// with box
template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, typename Interior, bool granularity, typename... Args>
Node* BaseTree<Point, DerivedTree, kBDO>::RebuildSingleTree(Node* T,
                                                            Args&&... args) {
    Points wx, wo;
    PrepareRebuild<Leaf, Interior, granularity>(T, wx, wo);
    static_assert(std::is_invocable_v<decltype(&DerivedTree::BuildRecursive),
                                      DerivedTree*, Slice, Slice, Args&&...>);
    return static_cast<DerivedTree*>(this)->BuildRecursive(
        parlay::make_slice(wx), parlay::make_slice(wo),
        std::forward<Args>(args)...);
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<SupportsForceParallel Interior, bool granularity>
inline bool BaseTree<Point, DerivedTree, kBDO>::ForceParallelRecursion(
    const Interior* TI) {
    return (granularity && TI->size > kSerialBuildCutoff) ||
           (!granularity && TI->ForceParallel());
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsBinaryNode Interior>
Node* BaseTree<Point, DerivedTree, kBDO>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq& pivots,
    parlay::sequence<Node*>& tree_nodes) {
    if (idx > kPivotNum) {
        assert(idx - kPivotNum - 1 < kBucketNum);
        return tree_nodes[idx - kPivotNum - 1];
    }
    Node *L, *R;
    L = BuildInnerTree<Interior>(idx << 1, pivots, tree_nodes);
    R = BuildInnerTree<Interior>(idx << 1 | 1, pivots, tree_nodes);
    return AllocInteriorNode<Interior>(L, R, pivots[idx],
                                       typename Interior::AT());
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, IsBinaryNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kBDO>::FlattenRec(Node* T, Range Out) {
    assert(T->size == Out.size());

    if (T->size == 0) return;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        for (int i = 0; i < TL->size; i++) {
            Out[i] = TL->pts[(!TL->is_dummy) * i];
        }
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

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kBDO>::FlattenRec(Node* T, Range Out) {
    if (T->size != Out.size()) {
        std::cout << "T->size: " << T->size << " Out.size(): " << Out.size()
                  << std::endl;
    }
    assert(T->size == Out.size());

    if (T->size == 0) return;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        for (int i = 0; i < TL->size; i++) {
            Out[i] = TL->pts[(!TL->is_dummy) * i];
        }
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
                TI->tree_nodes[i],
                Out.cut(start, start + TI->tree_nodes[i]->size));
        });
    } else {
        size_t start = 0;
        for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
            FlattenRec<Leaf, Interior, Range>(
                TI->tree_nodes[i],
                Out.cut(start, start + TI->tree_nodes[i]->size));
            start += TI->tree_nodes[i]->size;
        }
    }

    return;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kBDO>::PartialFlatten(Node* T, Range Out,
                                                        BucketType idx) {

    if (idx == 1) {
        assert(T->size == Out.size());
        FlattenRec<Leaf, Interior>(T, Out.cut(0, T->size));
        return;
    } else if (idx >= Interior::kRegions) {
        Node* ns =
            static_cast<Interior*>(T)->tree_nodes[idx - Interior::kRegions];
        assert(ns->size == Out.size());
        FlattenRec<Leaf, Interior>(ns, Out.cut(0, ns->size));
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    size_t l_size = TI->ReduceSums(idx << 1),
           r_size = TI->ReduceSums(idx << 1 | 1);
    assert(l_size + r_size == Out.size());
    parlay::par_do_if(
        ForceParallelRecursion<Interior, granularity>(
            static_cast<Interior*>(T)),
        [&]() {
            PartialFlatten<Leaf, Interior>(T, Out.cut(0, l_size), idx << 1);
        },
        [&]() {
            PartialFlatten<Leaf, Interior>(T, Out.cut(l_size, l_size + r_size),
                                           idx << 1 | 1);
        });
    return;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<Point, DerivedTree, kBDO>::ExpandMultiNode(
    const typename MN::ST& split, BucketType idx, BucketType deep,
    const parlay::sequence<Node*>& tree_nodes) {
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

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsBinaryNode BN, IsMultiNode MN>
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>
Node* BaseTree<Point, DerivedTree, kBDO>::Expand2Binary(Node* T) {
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
template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsBinaryNode Interior>
inline void BaseTree<Point, DerivedTree, kBDO>::UpdateInterior(Node* T, Node* L,
                                                               Node* R) {
    assert(!T->is_leaf);
    Interior* TI = static_cast<Interior*>(T);
    TI->ResetParallelFlag();
    TI->size = L->size + R->size;
    TI->left = L;
    TI->right = R;
    return;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsBinaryNode Interior>
inline void BaseTree<Point, DerivedTree, kBDO>::UpdateInterior(
    Node* T, const NodeBox& L, const NodeBox& R) {
    assert(!T->is_leaf);
    Interior* TI = static_cast<Interior*>(T);
    TI->ResetParallelFlag();
    TI->size = L.first->size + R.first->size;
    TI->left = L.first;
    TI->right = R.first;
    return;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<IsMultiNode Interior>
inline void BaseTree<Point, DerivedTree, kBDO>::UpdateInterior(
    Node* T, const typename Interior::NodeArr& new_nodes) {
    assert(!T->is_leaf);
    Interior* TI = static_cast<Interior*>(T);
    TI->ResetParallelFlag();
    TI->size = std::accumulate(
        new_nodes.begin(), new_nodes.end(), 0,
        [](size_t acc, Node* n) -> size_t { return acc + n->size; });
    TI->tree_nodes = new_nodes;
    return;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf>
Node* BaseTree<Point, DerivedTree, kBDO>::InsertPoints2Leaf(Node* T, Slice In) {
    Leaf* TL = static_cast<Leaf*>(T);
    if (TL->is_dummy) {
        T->size += In.size();
        return T;
    }

    if (TL->pts.size() == 0) {
        TL->pts = Points::uninitialized(kLeaveWrap);
    }
    for (int i = 0; i < In.size(); i++) {
        TL->pts[TL->size + i] = In[i];
    }
    TL->size += In.size();
    return T;
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, typename RT>
RT BaseTree<Point, DerivedTree, kBDO>::DeletePoints4Leaf(Node* T, Slice In) {
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
            return NodeBox(
                T, T->size ? Box(TL->pts[0], TL->pts[0]) : GetEmptyBox());
        } else {
            static_assert(std::same_as<RT, Node*> || std::same_as<RT, NodeBox>);
        }
    }

    auto it = TL->pts.begin(), end = TL->pts.begin() + TL->size;
    for (int i = 0; i < In.size(); i++) {
        it = std::ranges::find(TL->pts.begin(), end, In[i]);
        assert(it != end);
        std::ranges::iter_swap(it, --end);
    }

    assert(std::distance(TL->pts.begin(), end) == TL->size - In.size());
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
}  // namespace cpdd
