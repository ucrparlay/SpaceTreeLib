#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point, uint_fast8_t kBDO>
template<typename Leaf, typename Interior>
struct BaseTree<Point, kBDO>::InnerTree {
    // NOTE: helpers
    bool AssertSize(Node* T) const {
        if (T->is_leaf) {
            Leaf* TI = static_cast<Leaf*>(T);
            assert(T->size <= TI->pts.size() && T->size <= kLeaveWrap);
            return true;
        }
        Interior* TI = static_cast<Interior*>(T);
        assert(TI->size == TI->left->size + TI->right->size);
        return true;
    }

    void AssertSizeByIdx(BucketType idx) const {
        if (idx > kPivotNum || tags[idx].first->is_leaf) return;
        Interior* TI = static_cast<Interior*>(tags[idx].first);
        assert(TI->size == TI->left->size + TI->right->size);
        return;
    }

    inline BucketType GetNodeIdx(BucketType idx, Node* T) {
        if (tags[idx].first == T) return idx;
        if (idx > kPivotNum || tags[idx].first->is_leaf) return -1;
        auto pos = GetNodeIdx(idx << 1, T);
        if (pos != -1) return pos;
        return GetNodeIdx(idx << 1 | 1, T);
    }

    // NOTE: cores
    inline void ResetTagsNum() { tags_num = 0; }

    // NOTE: Each Node in the skeleton receives a tag
    // NOTE: A Leaf Node receives the tag < BUCKETNUM
    // NOTE: All internal Node has tag == BUCKETNUM
    void AssignNodeTag(Node* T, BucketType idx) {
        if (T->is_leaf || idx > kPivotNum) {
            assert(tags_num < kBucketNum);
            tags[idx] = NodeTag(T, tags_num);
            rev_tag[tags_num++] = idx;  // WARN: cannot remove
            return;
        }
        // INFO: BUCKET ID in [0, kBucketNum)
        tags[idx] = NodeTag(T, kBucketNum);
        Interior* TI = static_cast<Interior*>(T);
        if constexpr (IsBinaryNode<Interior>) {
            AssignNodeTag(TI->left, idx << 1);
            AssignNodeTag(TI->right, idx << 1 | 1);
        } else if constexpr (IsMultiNode<Interior>) {
            for (BucketType i = 0; i < Interior::kRegions; ++i) {
                AssignNodeTag(TI->tree_nodes[i], idx * Interior::kRegions + i);
            }
        } else {
            ;
        }
        return;
    }

    void ReduceSums(BucketType idx) const {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            assert(tags[idx].second < kBucketNum);
            sums_tree[idx] = sums[tags[idx].second];
            return;
        }
        ReduceSums(idx << 1);
        ReduceSums(idx << 1 | 1);
        sums_tree[idx] = sums_tree[idx << 1] + sums_tree[idx << 1 | 1];
        return;
    }

    static inline DimsType GetDepthByIndex(BucketType tag) {
        DimsType h = 0;
        while (tag > 1) {
            tag >>= 1;
            ++h;
        }
        return h;
    }

    // NOTE: a bucket/leaf has id kBucketNum+1
    // a node needs to be rebuilt has id kBucketNum+2
    // otherwise, it has id kBucketNum
    void PickTag(BucketType idx) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            tags[idx].second = kBucketNum + 1;
            rev_tag[tags_num++] = idx;
            return;
        }
        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        Interior* TI = static_cast<Interior*>(tags[idx].first);

        if (ImbalanceNode(TI->left->size + sums_tree[idx << 1],
                          TI->size + sums_tree[idx])) {
            tags[idx].second = kBucketNum + 2;
            rev_tag[tags_num++] = idx;
            return;
        }
        PickTag(idx << 1);
        PickTag(idx << 1 | 1);
        return;
    }

    void TagInbalanceNode() {
        ReduceSums(1);
        ResetTagsNum();
        PickTag(1);
        assert(AssertSize(tags[1].first));
        return;
    }

    void RetriveRebuildTreeIdx(BucketType idx,
                               parlay::sequence<BucketType>& re_idx,
                               size_t& tot_re_sz, BucketType& p) {
        if (tags[idx].second == kBucketNum + 3) {
            tot_re_sz += tags[idx].first->size;
            re_idx[p++] = idx;
        }
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            return;
        }
        RetriveRebuildTreeIdx(idx << 1, re_idx, tot_re_sz, p);
        RetriveRebuildTreeIdx(idx << 1 | 1, re_idx, tot_re_sz, p);
        return;
    }

    // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
    // the bucket Node whose ancestor has been rebuilt has tag kBucketNum+2
    // the bucket Node whose ancestor has not been ... has kBucketNum+1
    // otherwise, it's kBucketNum
    void MarkTomb(BucketType idx, BoxSeq& box_seq, Box bx, bool hasTomb) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            assert(tags[idx].second >= 0 && tags[idx].second < kBucketNum);
            if (!hasTomb) {  // NOTE: this subtree needs to be rebuilt in the
                             // future, therefore, ensure the granularity
                             // control by assign to aug_flag
                tags[idx].second = kBucketNum + 2;
                if (!tags[idx].first->is_leaf) {
                    // TODO: use api to assign parallel tag
                    assert(static_cast<Interior*>(tags[idx].first)
                               ->ForceParallel() == false);
                    static_cast<Interior*>(tags[idx].first)
                        ->SetParallelFlag(tags[idx].first->size >
                                          kSerialBuildCutoff);
                }
            } else {
                tags[idx].second = kBucketNum + 1;
            }
            box_seq[tags_num] = bx;
            rev_tag[tags_num++] = idx;
            return;
        }

        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        Interior* TI = static_cast<Interior*>(tags[idx].first);
        if (hasTomb && (ImbalanceNode(TI->left->size - sums_tree[idx << 1],
                                      TI->size - sums_tree[idx]) ||
                        (TI->size - sums_tree[idx] < kThinLeaveWrap))) {
            assert(hasTomb != 0);
            assert(TI->ForceParallel() == 0);
            tags[idx].second = kBucketNum + 3;
            hasTomb = false;
        }

        // NOTE: hasTomb == false => need to rebuild
        TI->SetParallelFlag(hasTomb ? false : TI->size > kSerialBuildCutoff);

        // TODO: rewrite for write efficient manner
        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        MarkTomb(idx << 1, box_seq, lbox, hasTomb);
        MarkTomb(idx << 1 | 1, box_seq, rbox, hasTomb);
        return;
    }

    void TagInbalanceNodeDeletion(BoxSeq& box_seq, Box bx, bool hasTomb) {
        ReduceSums(1);
        ResetTagsNum();
        MarkTomb(1, box_seq, bx, hasTomb);
        assert(AssertSize(tags[1].first));
        return;
    }

    InnerTree()
        requires IsBinaryNode<Interior>
        :
        tags_num(0),
        tags(NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1)),
        sums_tree(parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1)),
        rev_tag(Tag2Node::uninitialized(kBucketNum)) {}

    InnerTree()
        requires IsMultiNode<Interior>
        :
        tags_num(0),
        tags(NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1)),
        rev_tag(Tag2Node::uninitialized(kBucketNum)) {}

    // enum UpdateType { kPointer, kPointerBox };

    template<typename Base, typename Func>
    Base UpdateInnerTree(parlay::sequence<Base>& tree_nodes, Func&& func) {
        BucketType p = 0;
        return UpdateInnerTreeRecursive(1, tree_nodes, p,
                                        std::forward<Func>(func));
    }

    template<typename Base, typename Func>
    Base UpdateInnerTreeRecursive(BucketType idx,
                                  parlay::sequence<Base>& tree_nodes,
                                  BucketType& p, Func&& func) {
        if (this->tags[idx].second == kBucketNum + 1 ||
            this->tags[idx].second == kBucketNum + 2) {
            assert(this->rev_tag[p] == idx);
            return tree_nodes[p++];
        }

        // TODO:: perf
        Base left = UpdateInnerTreeRecursive(idx << 1, tree_nodes, p, func);
        Base right =
            UpdateInnerTreeRecursive(idx << 1 | 1, tree_nodes, p, func);

        static_assert(
            std::is_same_v<Base,
                           decltype(func(left, right, this->tags[idx], idx))>);
        return func(left, right, this->tags[idx], idx);
    }

    void Reset() {
        ResetTagsNum();
        tags = NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1);
        sums_tree = parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1);
        rev_tag = Tag2Node::uninitialized(kBucketNum);
    }

    // NOTE: variables
    NodeTagSeq tags;  //@ Assign each node a tag, aka skeleton
    parlay::sequence<BallsType> sums;
    mutable parlay::sequence<BallsType> sums_tree;
    mutable Tag2Node rev_tag;  //@ maps tag to the position in skeleton
    BucketType tags_num;
};

};  // namespace cpdd
