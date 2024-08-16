#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point, uint_fast8_t kBDO>
template<typename Leaf, typename Interior>
struct BaseTree<Point, kBDO>::InnerTree {
    //@ helpers
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

    //@ cores
    inline void ResetTagsNum() { tagsNum = 0; }

    // NOTE: Each Node in the skeleton receives a tag
    // NOTE: A Leaf Node receives the tag < BUCKETNUM
    // NOTE: All internal Node has tag == BUCKETNUM
    void AssignNodeTag(Node* T, BucketType idx) {
        if (T->is_leaf || idx > kPivotNum) {
            assert(tagsNum < kBucketNum);
            tags[idx] = NodeTag(T, tagsNum);
            rev_tag[tagsNum++] = idx;
            return;
        }
        // INFO: BUCKET ID in [0, kBucketNum)
        tags[idx] = NodeTag(T, kBucketNum);
        Interior* TI = static_cast<Interior*>(T);
        AssignNodeTag(TI->left, idx << 1);
        AssignNodeTag(TI->right, idx << 1 | 1);
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

    DimsType GetDepthByIndex(BucketType tag) {
        DimsType h = 0;
        while (tag > 1) {
            tag >>= 1;
            ++h;
        }
        return h;
    }

    void PickTag(BucketType idx) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            tags[idx].second = kBucketNum + 1;
            rev_tag[tagsNum++] = idx;
            return;
        }
        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        Interior* TI = static_cast<Interior*>(tags[idx].first);

        if (ImbalanceNode(TI->left->size + sums_tree[idx << 1],
                          TI->size + sums_tree[idx])) {
            tags[idx].second = kBucketNum + 2;
            rev_tag[tagsNum++] = idx;
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

    // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
    // the bucket Node whose ancestor has been rebuilt has tag kBucketNum+2
    // the bucket Node whose ancestor has not been ... has kBucketNum+1
    // otherwise, it's kBucketNum
    void MarkTomb(BucketType idx, BoxSeq& boxs, Box bx, bool hasTomb) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            assert(tags[idx].second >= 0 && tags[idx].second < kBucketNum);
            // tags[idx].second = (!hasTomb) ? kBucketNum + 2 : kBucketNum + 1;
            if (!hasTomb) {  // NOTE: this subtree needs to be rebuilt in the
                             // future, therefore ensure the granularity control
                             // by assign to aug_flag
                tags[idx].second = kBucketNum + 2;
                if (!tags[idx].first->is_leaf) {
                    assert(static_cast<Interior*>(tags[idx].first)->aug_flag =
                               false);
                    static_cast<Interior*>(tags[idx].first)->aug_flag =
                        tags[idx].first->size > kSerialBuildCutoff;
                }
            } else {
                tags[idx].second = kBucketNum + 1;
            }
            boxs[tagsNum] = bx;
            rev_tag[tagsNum++] = idx;
            return;
        }

        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        Interior* TI = static_cast<Interior*>(tags[idx].first);
        if (hasTomb && (inbalance_node(TI->left->size - sums_tree[idx << 1],
                                       TI->size - sums_tree[idx]) ||
                        (TI->size - sums_tree[idx] < kThinLeaveWrap))) {
            assert(hasTomb != 0);
            assert(TI->aug_flag == 0);
            tags[idx].second = kBucketNum + 3;
            hasTomb = false;
        }

        // NOTE: hasTomb == false => need to rebuild
        TI->aug_flag = hasTomb ? false : TI->size > kSerialBuildCutoff;

        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        MarkTomb(idx << 1, boxs, lbox, hasTomb);
        MarkTomb(idx << 1 | 1, boxs, rbox, hasTomb);
        return;
    }

    void TagInbalanceNodeDeletion(BoxSeq& boxs, Box bx, bool hasTomb) {
        ReduceSums(1);
        ResetTagsNum();
        MarkTomb(1, boxs, bx, hasTomb);
        assert(AssertSize(tags[1].first));
        return;
    }

    void Init() {
        ResetTagsNum();
        tags = NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1);
        sums_tree = parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1);
        rev_tag = TagNodes::uninitialized(kBucketNum);
    }

    //@ variables
    NodeTagSeq tags;
    parlay::sequence<BallsType> sums;
    mutable parlay::sequence<BallsType> sums_tree;
    mutable TagNodes rev_tag;
    BucketType tagsNum;
};

};  // namespace cpdd
