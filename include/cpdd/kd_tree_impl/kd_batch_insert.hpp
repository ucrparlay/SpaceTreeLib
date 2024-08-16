#pragma once

#include "../kd_tree.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchInsert(Slice A) {
    if (this->root_ == nullptr) {
        return Build_(A);  // BUG:check
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
    BucketType idx, const NodeTagSeq& tags, parlay::sequence<Node*>& treeNodes,
    BucketType& p, const TagNodes& rev_tag) {

    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        assert(rev_tag[p] == idx);
        return treeNodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != NULL);
    Node *L, *R;
    L = UpdateInnerTreeByTag(idx << 1, tags, treeNodes, p, rev_tag);
    R = UpdateInnerTreeByTag(idx << 1 | 1, tags, treeNodes, p, rev_tag);
    UpdateInterior(tags[idx].first, L, R);
    return tags[idx].first;
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::RebuildWithInsert(Node* T, Slice In,
                                                        const DimsType d) {
    DimsType curDim = this->split_rule_.FindRebuildDimension(d);
    Points wo = Points::uninitialized(T->size + In.size());
    Points wx = Points::uninitialized(T->size + In.size());
    parlay::parallel_for(0, In.size(), [&](size_t j) { wx[j] = In[j]; });
    BT::template FlattenRec<Leaf, Interior>(T, wx.cut(In.size(), wx.size()));
    BT::template DeleteTreeRecursive<Leaf, Interior>(T);
    return BuildRecursive(parlay::make_slice(wx), parlay::make_slice(wo),
                          curDim, BT::GetBox(parlay::make_slice(wx)));
}

//* return the updated Node
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
        UpdateInterior(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);
        return T;
    }

    // NOTE: assign each Node a tag
    typename BT::template InnerTree<Leaf, Interior> IT;
    IT.Init();
    assert(IT.rev_tag.size() == BT::kBucketNum);
    IT.AssignNodeTag(T, 1);
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);

    SeievePoints(In, Out, n, IT.tags, IT.sums, IT.tagsNum);

    IT.TagInbalanceNode();
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);
    auto treeNodes = parlay::sequence<Node*>::uninitialized(IT.tagsNum);

    parlay::parallel_for(
        0, IT.tagsNum,
        [&](size_t i) {
            size_t s = 0;
            for (int j = 0; j < i; j++) {
                s += IT.sums_tree[IT.rev_tag[j]];
            }

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
            if (IT.tags[IT.rev_tag[i]].second ==
                BT::kBucketNum + 1) {  // NOTE: continue sieve
                treeNodes[i] = BatchInsertRecursive(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
                    In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            } else {  // NOTE: launch rebuild subtree
                assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 2);
                assert(IT.tags[IT.rev_tag[i]].first->size +
                           IT.sums_tree[IT.rev_tag[i]] >=
                       0);

                treeNodes[i] = RebuildWithInsert(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            }
        },
        1);

    BucketType beatles = 0;
    return UpdateInnerTreeByTag(1, IT.tags, treeNodes, beatles, IT.rev_tag);
}
// NOTE: update the info of T by new children L and R
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
inline void KdTree<Point, SplitRule, kBDO>::UpdateInterior(Node* T, Node* L,
                                                           Node* R) {
    assert(!T->is_leaf);
    Interior* TI = static_cast<Interior*>(T);
    TI->size = L->size + R->size;
    TI->left = L;
    TI->right = R;
    return;
}

// NOTE: retrive the bucket tag of Point p from the skeleton tags
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::BucketType
KdTree<Point, SplitRule, kBDO>::RetriveTag(const Point& p,
                                           const NodeTagSeq& tags) {
    BucketType k = 1;
    Interior* TI;
    while (k <= BT::kPivotNum && (!tags[k].first->is_leaf)) {
        TI = static_cast<Interior*>(tags[k].first);
        k = Num::Lt(p.pnt[TI->split.second], TI->split.first) ? k << 1
                                                              : k << 1 | 1;
    }
    assert(tags[k].second < BT::kBucketNum);
    return tags[k].second;
}

// NOTE: seieve Points from range A to range B, using the skeleton tags. The
// sums is the number of elemenets within each bucket, the tagsNum is the total
// number of buckets in the skeleton
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::SeievePoints(
    Slice A, Slice B, const size_t n, const NodeTagSeq& tags,
    parlay::sequence<BallsType>& sums, const BucketType tagsNum) {
    size_t num_block = (n + BT::kBlockSize - 1) >> BT::kLog2Base;
    parlay::sequence<parlay::sequence<BallsType>> offset(
        num_block, parlay::sequence<BallsType>(tagsNum));
    assert(offset.size() == num_block && offset[0].size() == tagsNum &&
           offset[0][0] == 0);
    parlay::parallel_for(0, num_block, [&](size_t i) {
        for (size_t j = i << BT::kLog2Base;
             j < std::min((i + 1) << BT::kLog2Base, n); j++) {
            offset[i][std::move(RetriveTag(A[j], tags))]++;
        }
    });

    sums = parlay::sequence<BallsType>(tagsNum);
    for (size_t i = 0; i < num_block; i++) {
        auto t = offset[i];
        offset[i] = sums;
        for (int j = 0; j < tagsNum; j++) {
            sums[j] += t[j];
        }
    }

    parlay::parallel_for(0, num_block, [&](size_t i) {
        auto v = parlay::sequence<BallsType>::uninitialized(tagsNum);
        int tot = 0, s_offset = 0;
        for (int k = 0; k < tagsNum - 1; k++) {
            v[k] = tot + offset[i][k];
            tot += sums[k];
            s_offset += offset[i][k];
        }
        v[tagsNum - 1] = tot + ((i << BT::kLog2Base) - s_offset);
        for (size_t j = i << BT::kLog2Base;
             j < std::min((i + 1) << BT::kLog2Base, n); j++) {
            B[v[std::move(RetriveTag(A[j], tags))]++] = A[j];
        }
    });

    return;
}

// NOTE: traverse the skeleton tags recursively and update its children to new
// ones
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::UpdateInnerTree(
    BucketType idx, const NodeTagSeq& tags,
    parlay::sequence<NodeBox>& treeNodes, BucketType& p,
    const TagNodes& rev_tag) {
    if (tags[idx].second < BT::kBucketNum) {
        assert(rev_tag[p] == idx);
        return treeNodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != nullptr);
    auto [L, Lbox] = UpdateInnerTree(idx << 1, tags, treeNodes, p, rev_tag);
    auto [R, Rbox] = UpdateInnerTree(idx << 1 | 1, tags, treeNodes, p, rev_tag);
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

    UpdateInterior(T, L, R);

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
}

}  // namespace cpdd
