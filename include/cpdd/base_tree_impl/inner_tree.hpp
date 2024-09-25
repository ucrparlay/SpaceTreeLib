#pragma once

#include <utility>
#include "../base_tree.h"

namespace cpdd {
template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
template<typename Leaf, typename Interior>
struct BaseTree<Point, DerivedTree, kBDO>::InnerTree {
    using BT = BaseTree<Point, DerivedTree, kBDO>;
    InnerTree(BT& _BTRef) :
        BTRef(_BTRef),
        tags_num(0),
        tags(NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1)),
        sums_tree(parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1)),
        rev_tag(BucketSeq::uninitialized(kBucketNum)) {}

    // NOTE: helpers
    bool AssertSize(Node* T) const {
        if (T->is_leaf) {
            Leaf* TI = static_cast<Leaf*>(T);
            assert(T->size <= TI->pts.size() && T->size <= kLeaveWrap);
            return true;
        }
        Interior* TI = static_cast<Interior*>(T);
        if constexpr (IsBinaryNode<Interior>) {
            assert(TI->size == TI->left->size + TI->right->size);
        } else if constexpr (IsMultiNode<Interior>) {
            assert(TI->size ==
                   std::accumulate(
                       TI->tree_nodes.begin(), TI->tree_nodes.end(), 0,
                       [](size_t sum, Node* T) { return sum + T->size; }));
        } else {
            static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
        }
        return true;
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

    Box GetBoxByIdx(BucketType idx, const Box& box)
        requires IsMultiNode<Interior>
    {
        assert(GetDepthByIndex(19) == 4);
        assert(kDim = std::log2(Interior::kRegions));

        Box bx(box);
        BucketType h = GetDepthByIndex(idx);
        // idx -= (1 << h);
        // for (BucketType i = h, new_idx = 1; i > 0; i -= kDim) {
        //     auto TI = static_cast<Interior*>(tags[new_idx].first);
        //     BucketType local_id = 1;
        //     for (BucketType j = 0; j < kDim; j++) {
        //         bool bit_set = idx & (1 << (i - 1 - j));
        //         local_id = (local_id << 1) | bit_set;
        //         new_idx = (new_idx << 1) | bit_set;
        //     }
        //     TI->ModifyBoxById(local_id - Interior::kRegions, bx);
        // }

        for (BucketType i = h, new_idx = 1; i > 0; i -= kDim) {
            BucketType local_id = (idx >> (i - kDim)) & ((1 << kDim) - 1);
            auto TI = static_cast<Interior*>(tags[new_idx].first);
            TI->ModifyBoxById(local_id, bx);
            new_idx = new_idx << kDim | local_id;
        }
        return std::move(bx);
    }

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
            static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
        }
        return;
    }

    void ReduceSums(BucketType idx) const {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            assert(tags[idx].second < kBucketNum);
            sums_tree[idx] = sums[tags[idx].second];
            return;
        }

        if constexpr (IsBinaryNode<Interior>) {
            ReduceSums(idx << 1);
            ReduceSums(idx << 1 | 1);
            sums_tree[idx] = sums_tree[idx << 1] + sums_tree[idx << 1 | 1];
        } else if constexpr (IsMultiNode<Interior>) {
            sums_tree[idx] = 0;
            for (BucketType i = 0; i < Interior::kRegions; ++i) {
                ReduceSums(idx * Interior::kRegions + i);
                sums_tree[idx] += sums_tree[idx * Interior::kRegions + i];
            }
        }
        return;
    }

    static inline BucketType GetDepthByIndex(BucketType idx) {
        BucketType h = 0;
        while (idx > 1) {
            idx >>= 1;
            ++h;
        }
        return h;
    }

    // NOTE: a bucket/leaf has id kBucketNum+1
    // a node needs to be rebuilt has id kBucketNum+2
    // otherwise, it has id kBucketNum
    template<typename ViolateFunc>
        requires std::is_invocable_r_v<bool, ViolateFunc, BucketType>
    void PickTag(BucketType idx, ViolateFunc&& violate_func) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            tags[idx].second = kBucketNum + 1;
            rev_tag[tags_num++] = idx;
            return;
        }

        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        if (violate_func(idx)) {
            tags[idx].second = kBucketNum + 2;
            rev_tag[tags_num++] = idx;
            return;
        }

        if constexpr (IsBinaryNode<Interior>) {
            PickTag(idx << 1, violate_func);
            PickTag(idx << 1 | 1, violate_func);
        } else if constexpr (IsMultiNode<Interior>) {
            for (BucketType i = 0; i < Interior::kRegions; ++i) {
                PickTag(idx * Interior::kRegions + i, violate_func);
            }
        }
        return;
    }

    template<typename... Args>
    void TagInbalanceNode(Args&&... args) {
        ReduceSums(1);
        ResetTagsNum();
        PickTag(1, std::forward<Args>(args)...);
        assert(AssertSize(tags[1].first));
        return;
    }

    // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
    // the *bucket* Node whose ancestor has been rebuilt has tag kBucketNum+2
    // the *bucket* Node whose ancestor has not been ... has kBucketNum+1
    // otherwise, it's kBucketNum
    template<typename ViolateFunc>
        requires std::predicate<ViolateFunc, BucketType>
    void MarkTomb(BucketType idx, BucketType& re_num, size_t& tot_re_size,
                  BoxSeq& box_seq, const Box& box, bool has_tomb,
                  ViolateFunc&& violat_func) {
        if (idx > kPivotNum || tags[idx].first->is_leaf) {
            assert(tags[idx].second >= 0 && tags[idx].second < kBucketNum);
            if (!has_tomb) {  // NOTE: this subtree needs to be rebuilt in the
                              // future, therefore, ensure the granularity
                              // control by assign to aug_flag
                tags[idx].second = kBucketNum + 2;
                if (!tags[idx].first->is_leaf) {
                    auto TI = static_cast<Interior*>(tags[idx].first);
                    assert(!TI->GetParallelFlagIniStatus());
                    TI->SetParallelFlag(TI->size > kSerialBuildCutoff);
                    assert(TI->GetParallelFlagIniStatus());
                }
            } else {
                tags[idx].second = kBucketNum + 1;
            }

            box_seq[tags_num] = box;
            rev_tag[tags_num++] = idx;
            return;
        }

        // NOTE: no need to mark the internal nodes with tag kBucketNum
        assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
        Interior* TI = static_cast<Interior*>(tags[idx].first);
        if (has_tomb && violat_func(idx)) {
            assert(has_tomb != 0);
            assert(!TI->GetParallelFlagIniStatus());
            tags[idx].second = kBucketNum + 3;
            has_tomb = false;
            TI->SetParallelFlag(TI->size > kSerialBuildCutoff);
            re_num++;
            tot_re_size += TI->size;
        }

        if constexpr (IsBinaryNode<Interior>) {
            // TODO: rewrite to reduce work
            Box lbox(box), rbox(box);
            lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
            rbox.first.pnt[TI->split.second] = TI->split.first;
            MarkTomb(idx << 1, re_num, tot_re_size, box_seq, lbox, has_tomb,
                     violat_func);
            MarkTomb(idx << 1 | 1, re_num, tot_re_size, box_seq, rbox, has_tomb,
                     violat_func);
        } else if constexpr (IsMultiNode<Interior>) {
            BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));
            for (BucketType i = 0; i < Interior::kRegions; ++i) {
                MarkTomb(idx * Interior::kRegions + i, re_num, tot_re_size,
                         box_seq, new_box[i], has_tomb, violat_func);
            }
        }
        return;
    }

    // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
    // the *bucket* Node whose ancestor has been rebuilt has tag kBucketNum+2
    // the *bucket* Node whose ancestor has not been ... has kBucketNum+1
    // otherwise, it's kBucketNum

    template<typename... Args>
    auto TagInbalanceNodeDeletion(Args&&... args) {
        ReduceSums(1);
        ResetTagsNum();
        BucketType re_num = 0;
        size_t tot_re_size = 0;
        MarkTomb(1, re_num, tot_re_size, std::forward<Args>(args)...);
        return std::make_pair(re_num, tot_re_size);
    }

    enum UpdateType { kUpdatePointer, kTagRebuildNode, kPostRebuild };

    // NOTE: update the skeleton based on the @UpdateType
    template<UpdateType kUT, typename Base, typename... Args>
        requires IsPointer<Base> ||
                 (IsPair<Base> && IsPointer<typename Base::first_type> &&
                  IsBox<typename Base::second_type, Point>)
    Base UpdateInnerTree(const parlay::sequence<Base>& tree_nodes,
                         Args&&... args) {
        BucketType p = 0;
        if constexpr (kUT ==
                      kUpdatePointer) {  // NOTE: update the inner tree nodes
            return UpdateInnerTreeRecursive<kUT>(1, tree_nodes, p, [&]() {});
        } else if constexpr (kUT ==
                             kTagRebuildNode) {  // NOTE: tag the node that
                                                 // needs to be reconstruct
            auto func = [&](const auto&... params) -> void {
                if constexpr (IsBinaryNode<Interior>) {  // needs to save the
                                                         // box
                    const auto& new_box =
                        std::get<0>(std::forward_as_tuple(params...));
                    const auto& idx =
                        std::get<1>(std::forward_as_tuple(params...));
                    rev_tag[tags_num] = idx;
                    FindVar<BoxSeq>(std::forward<Args>(args)...)[tags_num++] =
                        new_box;
                } else if constexpr (IsMultiNode<Interior>) {
                    const auto& idx =
                        std::get<0>(std::forward_as_tuple(params...));
                    rev_tag[tags_num++] = idx;
                }
            };

            this->ResetTagsNum();
            return UpdateInnerTreeRecursive<kUT>(1, tree_nodes, p, func);
        } else if constexpr (kUT ==
                             kPostRebuild) {  // NOTE: avoid touch the node
                                              // that has been deleted
            bool under_rebuild_tree = false;
            return UpdateInnerTreeRecursive<kUT>(
                1, tree_nodes, p, [&](bool op) -> bool {
                    return op == 0 ? (under_rebuild_tree = !under_rebuild_tree)
                                   : under_rebuild_tree;
                });
        } else {
            static_assert(0);
        }
    }

    // NOTE: udpate inner tree for binary nodes
    template<UpdateType kUT, typename Base, typename Func>
        requires IsBinaryNode<Interior>
    Base UpdateInnerTreeRecursive(BucketType idx,
                                  const parlay::sequence<Base>& tree_nodes,
                                  BucketType& p, Func&& func) {
        // WARN: needs to ensure this success for both insert and delete
        if (this->tags[idx].second == kBucketNum + 1 ||
            this->tags[idx].second == kBucketNum + 2) {
            // assert(this->rev_tag[p] == idx);
            return tree_nodes[p++];
        }

        if constexpr (kUT == kPostRebuild) {
            if (this->tags[idx].second == kBucketNum + 3) {
                func(0);  // close the under_rebuild_tree flag
                assert(func(1) == true);
            }
        }

        const Base& left =
            UpdateInnerTreeRecursive<kUT>(idx << 1, tree_nodes, p, func);
        const Base& right =
            UpdateInnerTreeRecursive<kUT>(idx << 1 | 1, tree_nodes, p, func);

        if constexpr (kUT ==
                      kUpdatePointer) {  // NOTE: only update the pointers
            UpdateInterior<Interior>(this->tags[idx].first, left, right);
            if constexpr (IsPointer<Base>) {
                return this->tags[idx].first;
            } else {  // WARN: if only update pointer, then avoid update box
                return NodeBox(this->tags[idx].first, Box());
            }
        } else if constexpr (kUT == kTagRebuildNode) {  // NOTE: retag and save
                                                        // box for rebuild
            UpdateInterior<Interior>(this->tags[idx].first, left, right);
            auto new_box = BT::GetBox(left.second, right.second);
            if (this->tags[idx].second == BT::kBucketNum + 3) {
                func(new_box, idx);
            }
            return NodeBox(this->tags[idx].first, std::move(new_box));
        } else if constexpr (kUT ==
                             kPostRebuild) {  // NOTE: avoid update pointers
                                              // for deleted trees
            if (!func(1)) {  // query whether under the rebuild_tree
                UpdateInterior<Interior>(this->tags[idx].first, left, right);
                return NodeBox(this->tags[idx].first,
                               Box());  // box has been computed before
            } else if (this->tags[idx].second == kBucketNum + 3) {  // back
                func(0);  // disable the under_rebuild_tree flag
                assert(func(1) == false);
                return NodeBox(this->tags[idx].first, Box());
            } else {  // the tree has been deleted
                return NodeBox(nullptr, Box());
            }
        } else {
            static_assert(0);
        }
    }

    template<UpdateType kUT, typename Base, typename Func>
        requires IsMultiNode<Interior>
    Base UpdateInnerTreeRecursive(BucketType idx,
                                  const parlay::sequence<Base>& tree_nodes,
                                  BucketType& p, Func&& func) {
        if (tags[idx].second == BT::kBucketNum + 1 ||
            tags[idx].second == BT::kBucketNum + 2) {
            return tree_nodes[p++];
        }

        if constexpr (kUT == kPostRebuild) {
            if (this->tags[idx].second == kBucketNum + 3) {
                func(0);  // close the under_rebuild_tree flag
                assert(func(1) == true);
            }
        }

        typename Interior::OrthNodeArr new_nodes;
        for (BucketType i = 0; i < Interior::kRegions; ++i) {
            new_nodes[i] = UpdateInnerTreeRecursive<kUT>(
                idx * Interior::kRegions + i, tree_nodes, p, func);
        }

        if constexpr (kUT == kUpdatePointer) {
            BT::template UpdateInterior<Interior>(tags[idx].first, new_nodes);
            return this->tags[idx].first;
        } else if constexpr (kUT == kTagRebuildNode) {
            UpdateInterior<Interior>(this->tags[idx].first, new_nodes);
            if (this->tags[idx].second == BT::kBucketNum + 3) {
                func(idx);
            }
            return this->tags[idx].first;
        } else if constexpr (kUT ==
                             kPostRebuild) {  // NOTE: avoid update pointers
                                              // for deleted trees
            if (!func(1)) {                   // not under rebuild tree
                UpdateInterior<Interior>(this->tags[idx].first, new_nodes);
                return this->tags[idx].first;  // box has been computed before
            } else if (this->tags[idx].second == kBucketNum + 3) {  // back
                func(0);  // disable the under_rebuild_tree flag
                assert(func(1) == false);
                return this->tags[idx].first;
            } else {  // the tree has been deleted
                return nullptr;
            }
        } else {
            static_assert(false);
        }
    }

    void Reset() {
        ResetTagsNum();
        tags = NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1);
        sums_tree = parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1);
        rev_tag = BucketSeq::uninitialized(kBucketNum);
    }

    // NOTE: variables
    NodeTagSeq tags;  //@ Assign each node a tag, aka skeleton
    parlay::sequence<BallsType> sums;
    mutable parlay::sequence<BallsType> sums_tree;
    mutable BucketSeq rev_tag;  //@ maps tag to the position in skeleton
    BucketType tags_num;
    BT& BTRef;
};
};  // namespace cpdd
