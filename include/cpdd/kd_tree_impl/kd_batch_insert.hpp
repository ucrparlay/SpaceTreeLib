#pragma once

#include "../kd_tree.h"
#include "parlay/slice.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchInsert(Slice A) {
    if (this->root_ == nullptr) {  // TODO: may check using explicity tag
        return Build_(A);
    }

    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    this->tree_box_ = BT::GetBox(this->tree_box_, BT::GetBox(A));
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    this->root_ = BatchInsertRecursive(T, A, B.cut(0, A.size()), d);
    assert(this->root_ != NULL);
    return;
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::UpdateInnerTreeByTag(
    BucketType idx, const NodeTagSeq& tags, parlay::sequence<Node*>& tree_nodes,
    BucketType& p, const TagNodes& rev_tag) {

    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        assert(rev_tag[p] == idx);
        return tree_nodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != NULL);
    Node *L, *R;
    L = UpdateInnerTreeByTag(idx << 1, tags, tree_nodes, p, rev_tag);
    R = UpdateInnerTreeByTag(idx << 1 | 1, tags, tree_nodes, p, rev_tag);
    BT::template UpdateInterior<Interior>(tags[idx].first, L, R);
    return tags[idx].first;
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::RebuildWithInsert(Node* T, Slice In,
                                                        const DimsType d) {
    DimsType curDim = this->split_rule_.FindRebuildDimension(d);
    Points wx, wo;
    BT::template PrepareRebuild<Leaf, Interior>(T, In, wx, wo);
    return BuildRecursive(parlay::make_slice(wx), parlay::make_slice(wo),
                          curDim, BT::GetBox(parlay::make_slice(wx)));
}

// NOTE: return the updated Node
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::BatchInsertRecursive(Node* T, Slice In,
                                                           Slice Out,
                                                           DimsType d) {
    size_t n = In.size();

    if (n == 0) return T;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        if (!TL->is_dummy && n + TL->size <= BT::kLeaveWrap) {
            assert(T->size == TL->size);
            if (TL->pts.size() == 0) {
                TL->pts = Points::uninitialized(BT::kLeaveWrap);
            }
            for (int i = 0; i < n; i++) {
                TL->pts[TL->size + i] = In[i];
            }
            TL->size += n;
            return T;
        } else {
            return RebuildWithInsert(T, In, d);
        }
    }

    if (n <= BT::kSerialBuildCutoff) {
        Interior* TI = static_cast<Interior*>(T);
        auto _2ndGroup = std::ranges::partition(In, [&](const Point& p) {
            return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });

        // NOTE: rebuild
        if (BT::ImbalanceNode(TI->left->size + _2ndGroup.begin() - In.begin(),
                              TI->size + n)) {
            return RebuildWithInsert(T, In, d);
        }

        // NOTE: continue
        Node *L, *R;
        d = (d + 1) % BT::kDim;
        L = BatchInsertRecursive(TI->left,
                                 In.cut(0, _2ndGroup.begin() - In.begin()),
                                 Out.cut(0, _2ndGroup.begin() - In.begin()), d);
        R = BatchInsertRecursive(TI->right,
                                 In.cut(_2ndGroup.begin() - In.begin(), n),
                                 Out.cut(_2ndGroup.begin() - In.begin(), n), d);
        BT::template UpdateInterior<Interior>(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);
        return T;
    }

    // NOTE: assign each Node a tag
    typename BT::template InnerTree<Leaf, Interior> IT;
    // IT.Init();
    assert(IT.rev_tag.size() == BT::kBucketNum);
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);

    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    IT.TagInbalanceNode();
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);

    parlay::parallel_for(
        0, IT.tags_num,
        [&](size_t i) {
            size_t s = 0;
            for (int j = 0; j < i; j++) {
                s += IT.sums_tree[IT.rev_tag[j]];
            }

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
            if (IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1) {
                // NOTE: continue sieve
                tree_nodes[i] = BatchInsertRecursive(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
                    In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            } else {  // NOTE: launch rebuild subtree
                assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 2);
                assert(IT.tags[IT.rev_tag[i]].first->size +
                           IT.sums_tree[IT.rev_tag[i]] >=
                       0);

                tree_nodes[i] = RebuildWithInsert(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            }
        },
        1);

    BucketType beatles = 0;
    return UpdateInnerTreeByTag(1, IT.tags, tree_nodes, beatles, IT.rev_tag);
}

// NOTE: traverse the skeleton tags and update its children to new ones
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::UpdateInnerTree(
    BucketType idx, const NodeTagSeq& tags,
    parlay::sequence<NodeBox>& tree_nodes, BucketType& p,
    const TagNodes& rev_tag) {
    if (tags[idx].second < BT::kBucketNum) {
        assert(rev_tag[p] == idx);
        return tree_nodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != nullptr);
    auto [L, Lbox] = UpdateInnerTree(idx << 1, tags, tree_nodes, p, rev_tag);
    auto [R, Rbox] =
        UpdateInnerTree(idx << 1 | 1, tags, tree_nodes, p, rev_tag);
    UpdateInterior(tags[idx].first, L, R);
    return NodeBox(tags[idx].first, BT::GetBox(Lbox, Rbox));
}

// NOTE: flatten a tree then rebuild upon it
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::RebuildSingleTree(Node* T, const DimsType d,
                                                  const bool granularity) {
    Points wo = Points::uninitialized(T->size);
    Points wx = Points::uninitialized(T->size);
    uint_fast8_t curDim = pick_rebuild_dim(T, d, BT::kDim);
    // flatten(T, wx.cut(0, T->size), granularity);
    // delete_tree_recursive(T, granularity);
    BT::template FlattenRec<Leaf, Interior>(T, wx.cut(0, T->size));
    BT::template DeleteTreeRecursive<Leaf, Interior>(T);
    Box bx = BT::GetBox(parlay::make_slice(wx));
    Node* o = build_recursive(parlay::make_slice(wx), parlay::make_slice(wo),
                              curDim, BT::kDim, bx);
    return NodeBox(std::move(o), std::move(bx));
}

// NOTE: traverse the tree in parallel and rebuild the imbalanced subtree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::RebuildTreeRecursive(Node* T, DimsType d,
                                                     const bool granularity) {
    if (T->is_leaf) {
        return NodeBox(T, BT::GetBox(T));
    }

    Interior* TI = static_cast<Interior*>(T);
    if (BT::ImBalanceNode(TI->left->size, TI->size)) {
        return RebuildSingleTree(T, d, BT::kDim, granularity);
    }

    Node *L, *R;
    Box Lbox, Rbox;
    d = (d + 1) % BT::kDim;
    parlay::par_do_if(
        // NOTE: if granularity is disabled, always traverse the tree in
        // parallel
        (granularity && T->size > BT::SerialBuildCutoff) ||
            (!granularity && TI->aug_flag),
        [&] {
            std::tie(L, Lbox) =
                RebuildTreeRecursive(TI->left, d, BT::kDim, granularity);
        },
        [&] {
            std::tie(R, Rbox) =
                RebuildTreeRecursive(TI->right, d, BT::kDim, granularity);
        });

    BT::template UpdateInterior<Interior>(T, L, R);

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
}

}  // namespace cpdd
