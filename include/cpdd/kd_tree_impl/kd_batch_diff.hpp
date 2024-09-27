#pragma once
#include "../kd_tree.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::BatchDiff(Range&& In) {
    static_assert(parlay::is_random_access_range_v<Range>);
    static_assert(parlay::is_less_than_comparable_v<
                  parlay::range_reference_type_t<Range>>);
    static_assert(
        std::is_constructible_v<parlay::range_value_type_t<Range>,
                                parlay::range_reference_type_t<Range>>);

    Slice A = parlay::make_slice(In);
    BatchDiff_(A);
    return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDiff_(Slice A) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    Box box = this->tree_box_;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    std::tie(T, this->tree_box_) =
        BatchDiffRecursive(T, box, A, parlay::make_slice(B), d);

    std::tie(this->root_, box) = RebuildTreeRecursive(T, d, false);
    // assert(box == this->tree_box_);

    return;
}

// NOTE: traverse the skeleton tags and update its children to new ones
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::UpdateInnerTreePointerBox(
    BucketType idx, const NodeTagSeq& tags, NodeBoxSeq& tree_nodes,
    BucketType& p) {
    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        return tree_nodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != nullptr);

    auto& [L, Lbox] = UpdateInnerTreePointerBox(idx << 1, tags, tree_nodes, p);
    auto& [R, Rbox] =
        UpdateInnerTreePointerBox(idx << 1 | 1, tags, tree_nodes, p);
    BT::template UpdateInterior<Interior>(tags[idx].first, L, R);
    return NodeBox(tags[idx].first, BT::GetBox(Lbox, Rbox));
}

// NOTE: traverse the tree in parallel and rebuild the imbalanced subtree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::RebuildTreeRecursive(Node* T, DimsType d,
                                                     const bool granularity) {
    if (T->is_leaf) {
        return NodeBox(T, BT::template GetBox<Leaf, Interior>(T));
    }

    Interior* TI = static_cast<Interior*>(T);
    if (BT::ImbalanceNode(TI->left->size, TI->size)) {
        // WARN: this disables the parallelism in default
        // TODO add box support
        return BT::template RebuildSingleTree<Leaf, Interior, false>(T, d);
    }

    Node *L, *R;
    Box Lbox, Rbox;
    d = (d + 1) % BT::kDim;
    parlay::par_do_if(
        // NOTE: if granularity is disabled, always traverse the tree in
        // parallel
        (granularity && T->size > BT::kSerialBuildCutoff) ||
            (!granularity && TI->aug),
        [&] {
            std::tie(L, Lbox) = RebuildTreeRecursive(TI->left, d, granularity);
        },
        [&] {
            std::tie(R, Rbox) = RebuildTreeRecursive(TI->right, d, granularity);
        });

    BT::template UpdateInterior<Interior>(T, L, R);

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
}

// NOTE: only sieve the Points, without rebuilding the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::BatchDiffRecursive(
    Node* T, const typename KdTree<Point, SplitRule, kBDO>::Box& box, Slice In,
    Slice Out, DimsType d) {
    size_t n = In.size();

    if (n == 0) return NodeBox(T, box);

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);

        if (TL->is_dummy) {  // NOTE: need to check whether all Points are in
                             // the Leaf
            // TODO: slow when In.size() is large
            for (size_t i = 0; TL->size && i < In.size(); i++) {
                if (TL->pts[0] == In[i]) {
                    TL->size--;
                }
            }
            assert(TL->size >= 0);
            return NodeBox(
                T, TL->size ? Box(TL->pts[0], TL->pts[0]) : BT::GetEmptyBox());
        }

        auto it = TL->pts.begin(), end = TL->pts.begin() + TL->size;
        for (int i = 0; TL->size && i < In.size(); i++) {
            it = std::ranges::find(TL->pts.begin(), end, In[i]);
            if (it != end) {  // NOTE: find a Point
                std::ranges::iter_swap(it, --end);
                TL->size--;
            }
        }
        return NodeBox(T, BT::GetBox(TL->pts.cut(0, TL->size)));
    }

    if (In.size() <= BT::kSerialBuildCutoff) {
        // if (In.size()) {
        Interior* TI = static_cast<Interior*>(T);
        PointsIter split_iter =
            std::ranges::partition(In, [&](const Point& p) {
                return Num::Lt(p.pnt[TI->split.second], TI->split.first);
            }).begin();

        DimsType nextDim = (d + 1) % BT::kDim;

        BoxCut box_cut(box, TI->split, true);
        auto& [L, Lbox] =
            BatchDiffRecursive(TI->left, box_cut.GetFirstBoxCut(),
                               In.cut(0, split_iter - In.begin()),
                               Out.cut(0, split_iter - In.begin()), nextDim);
        auto& [R, Rbox] =
            BatchDiffRecursive(TI->right, box_cut.GetSecondBoxCut(),
                               In.cut(split_iter - In.begin(), n),
                               Out.cut(split_iter - In.begin(), n), nextDim);

        BT::template UpdateInterior<Interior>(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        return NodeBox(T, BT::GetBox(Lbox, Rbox));
    }

    InnerTree IT(*this);
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
    auto boxs = parlay::sequence<Box>::uninitialized(IT.tags_num);

    // NOTE: never set tomb, this equivalent to only calcualte the bounding box,
    // BUG: cannot direct pass false here,
    // TODO: remove bounding boxes
    IT.TagInbalanceNodeDeletion(boxs, box, false);

    parlay::parallel_for(
        0, IT.tags_num,
        // NOTE: i is the index of the tags
        [&](size_t i) {
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
            tree_nodes[i] =
                BatchDiffRecursive(IT.tags[IT.rev_tag[i]].first, boxs[i],
                                   Out.cut(start, start + IT.sums[i]),
                                   In.cut(start, start + IT.sums[i]), nextDim);
        },
        1);

    BucketType beatles = 0;
    return UpdateInnerTreePointerBox(1, IT.tags, tree_nodes, beatles);
}
}  // namespace cpdd
