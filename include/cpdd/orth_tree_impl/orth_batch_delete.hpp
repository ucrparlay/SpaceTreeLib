#pragma once

#include <algorithm>
#include "../orth_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

// NOTE: default batch delete
template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
template<typename Range>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchDelete(Range&& In) {
    static_assert(parlay::is_random_access_range_v<Range>);
    static_assert(parlay::is_less_than_comparable_v<
                  parlay::range_reference_type_t<Range>>);
    static_assert(
        std::is_constructible_v<parlay::range_value_type_t<Range>,
                                parlay::range_reference_type_t<Range>>);

    Slice A = parlay::make_slice(In);
    BatchDelete_(A);
    return;
}

// NOTE: assume all Points are fully covered in the tree
template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchDelete_(Slice A) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    Box bx = this->tree_box_;
    std::tie(this->root_, this->tree_box_) =
        BatchDeleteRecursive(T, bx, A, parlay::make_slice(B), 1);
    return;
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::BatchDeleteRecursive(
    Node* T, Slice In, Slice Out, bool hasTomb) {
    size_t n = In.size();

    if (n == 0) {
        return T;
    }

    // INFO: may can be used to accelerate the whole deletion process
    if (n == T->size) {
        if (hasTomb) {
            BT::template DeleteTreeRecursive<Leaf, Interior>(T);
            return AllocEmptyLeafNode<Slice, Leaf>();
        }
        auto TI = static_cast<Interior*>(T);
        TI->SetParallelFlag(T->size > BT::kSerialBuildCutoff);
        TI->size = 0;
        return T;
    }

    if (T->is_leaf) {
        return BT::template DeletePoints4Leaf<Leaf, Node*>(T, In);
    }

    if (In.size() <= BT::kSerialBuildCutoff) {
        parlay::sequence<BallsType> sums(kNodeRegions, 0);
        SerialSplitSkeleton(T, In, 0, 1, sums);
        assert(std::accumulate(sums.begin(), sums.end(), 0) == n);

        // TODO: rebuild the tree if the current tree in the bb fallen behind
        // leave wrap
        Nodes new_nodes;
        size_t start = 0;
        for (DimsType i = 0; i < kNodeRegions; ++i) {
            new_nodes[i] =
                BatchDeleteRecursive(static_cast<Interior*>(T)->tree_nodes[i],
                                     In.cut(start, start + sums[i]),
                                     Out.cut(start, start + sums[i]));
            start += sums[i];
        }
        BT::template UpdateInterior<Interior>(T, new_nodes);
        assert(T->is_leaf == false);
        return T;
    }

    InnerTree IT(*this);
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
    auto box_seq = parlay::sequence<Box>::uninitialized(IT.tags_num);

    auto [re_num, tot_re_size] =
        IT.TagInbalanceNodeDeletion(box_seq, bx, hasTomb);

    // TODO: no need to compute, compute uing std::accumulate
    // IT.RetriveRebuildTreeIdx(1, re_idx, tot_re_size, re_num);
    assert(re_num <= IT.tags_num);

    parlay::parallel_for(
        0, IT.tags_num,
        [&](size_t i) {
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
            assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
            assert(BT::WithinBox(BT::GetBox(Out.cut(start, start + IT.sums[i])),
                                 BT::template GetBox<Leaf, Interior>(
                                     IT.tags[IT.rev_tag[i]].first)));

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;

            tree_nodes[i] = BatchDeleteRecursive(
                IT.tags[IT.rev_tag[i]].first, box_seq[i],
                Out.cut(start, start + IT.sums[i]),
                In.cut(start, start + IT.sums[i]), nextDim,
                IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
        },
        1);

    // NOTE: handling of rebuild
    // TODO: maybe we can use the tot_re_size/total_rebuild_num
    // > SERIAL_BUILD_CUTOFF to judge whether to rebuild the tree in parallel
    // WARN: the rebuild node is on top
    if (tot_re_size > BT::kSerialBuildCutoff) {  // NOTE: parallel rebuild

        // NOTE: get new box for skeleton root and rebuild nodes
        auto re_idx = BucketSeq::uninitialized(IT.tags_num);
        BucketType id = 0;  // PARA: the number of boxs needs to rebuild
        const Box new_box =
            std::get<1>(IT.template UpdateInnerTree<InnerTree::kSaveBox>(
                tree_nodes,
                [&](const Box& new_box, const BucketType idx) -> void {
                    re_idx[id] = idx;
                    box_seq[id++] = new_box;
                }));
        assert(id == re_num);

        parlay::parallel_for(0, re_num, [&](size_t i) {
            assert(IT.tags[re_idx[i]].second == BT::kBucketNum + 3);

            if (IT.tags[re_idx[i]].first->size == 0) {  // NOTE: empty
                BT::template DeleteTreeRecursive<Leaf, Interior, false>(
                    IT.tags[re_idx[i]].first);
                IT.tags[re_idx[i]].first = AllocEmptyLeafNode<Slice, Leaf>();
            } else {  // NOTE: rebuild
                auto next_dim = (d + IT.GetDepthByIndex(re_idx[id])) % BT::kDim;
                IT.tags[re_idx[i]].first = std::get<0>(
                    BT::template RebuildSingleTree<Leaf, Interior, false>(
                        IT.tags[re_idx[i]].first, d, box_seq[i]));
            }
        });
        bool under_rebuild_tree = false;
        const auto new_root =
            std::get<0>(IT.template UpdateInnerTree<InnerTree::kReturnRebuild>(
                tree_nodes, [&](bool op) -> bool {
                    return op == 0 ? (under_rebuild_tree = !under_rebuild_tree)
                                   : under_rebuild_tree;
                }));
        return NodeBox(new_root, new_box);
    } else {  // NOTE: update tree in serial
        return IT.template UpdateInnerTree<InnerTree::kDelete>(
            tree_nodes, [&](BucketType idx) -> DimsType {
                return (d + IT.GetDepthByIndex(idx)) % BT::kDim;
            });
    }
}

}  // namespace cpdd
