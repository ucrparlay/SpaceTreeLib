#pragma once

#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

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
// NOTE: default batch delete
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Range&& In) {
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
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete_(Slice A) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    Box bx = this->tree_box_;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    std::tie(this->root_, this->tree_box_) =
        BatchDeleteRecursive(T, bx, A, parlay::make_slice(B), d, 1);
    return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Slice A, PartialCoverTag) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    Box box = this->tree_box_;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    // NOTE: first sieve the Points
    std::tie(T, this->tree_box_) = BatchDeleteRecursive(
        T, box, A, parlay::make_slice(B), d, PartialCoverTag());
    // NOTE: then rebuild the tree with full parallelsim
    std::tie(this->root_, box) = RebuildTreeRecursive(T, d, false);
    assert(box == this->tree_box_);

    return;
}

// NOTE: the Node which needs to be rebuilt has tag BT::kBucketNum+3
// NOTE: the bucket Node whose ancestor has been rebuilt has tag
// BT::kBucketNum+2
// NOTE: the bucket Node whose ancestor has not been ... has
// BT::kBucketNum+1
// NOTE: otherwise, it's BT::kBucketNum
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::DeleteInnerTree(
    BucketType idx, const NodeTagSeq& tags,
    parlay::sequence<NodeBox>& tree_nodes, BucketType& p,
    const TagNodes& rev_tag, DimsType d) {
    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        assert(rev_tag[p] == idx);
        assert(tags[idx].second == BT::kBucketNum + 1 ||
               tags[idx].first->size > BT::kSerialBuildCutoff ==
                   static_cast<Interior*>(tags[idx].first)->aug);
        return tree_nodes[p++];  // WARN: this blocks the parallelsim as it will
                                 // rebuild the tree one-by-one
    }

    auto [L, Lbox] = DeleteInnerTree(idx << 1, tags, tree_nodes, p, rev_tag,
                                     (d + 1) % BT::kDim);
    auto [R, Rbox] = DeleteInnerTree(idx << 1 | 1, tags, tree_nodes, p, rev_tag,
                                     (d + 1) % BT::kDim);

    assert(tags[idx].first->size > BT::kSerialBuildCutoff ==
           static_cast<Interior*>(tags[idx].first)->aug);
    BT::template UpdateInterior<Interior>(tags[idx].first, L, R);

    if (tags[idx].second == BT::kBucketNum + 3) {  // NOTE: launch rebuild
        Interior const* TI = static_cast<Interior*>(tags[idx].first);
        assert(BT::ImbalanceNode(TI->left->size, TI->size) ||
               TI->size < BT::kThinLeaveWrap);

        if (tags[idx].first->size == 0) {  // NOTE: special judge for empty tree
            BT::template DeleteTreeRecursive<Leaf, Interior, false>(
                tags[idx].first);
            return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(),
                           BT::GetEmptyBox());
        }

        return BT::template RebuildSingleTree<Leaf, Interior, false>(
            tags[idx].first, d);
    }

    return NodeBox(tags[idx].first, BT::GetBox(Lbox, Rbox));
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
// WARN: the param d can be only used when rotate cutting is applied
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::BatchDeleteRecursive(
    Node* T, const typename KdTree<Point, SplitRule, kBDO>::Box& bx, Slice In,
    Slice Out, DimsType d, bool hasTomb) {
    size_t n = In.size();

    if (n == 0) {
        assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), bx));
        return NodeBox(T, bx);
    }

    // INFO: may can be used to accelerate the whole deletion process
    if (n == T->size) {
        if (hasTomb) {
            BT::template DeleteTreeRecursive<Leaf, Interior>(T);
            return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(),
                           BT::GetEmptyBox());
        }
        T->size = 0;
        return NodeBox(T, BT::GetEmptyBox());
    }

    if (T->is_leaf) {
        assert(T->size >= In.size());
        Leaf* TL = static_cast<Leaf*>(T);

        if (TL->is_dummy) {
            assert(T->is_leaf);
            assert(In.size() <=
                   T->size);  // WARN: cannot delete more Points then there are
            T->size -= In.size();  // WARN: this assumes that In\in T
            return NodeBox(T, Box(TL->pts[0], TL->pts[0]));
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
        return NodeBox(T, BT::GetBox(TL->pts.cut(0, TL->size)));
    }

    if (In.size() <= BT::kSerialBuildCutoff) {
        Interior* TI = static_cast<Interior*>(T);
        PointsIter split_iter =
            std::ranges::partition(In, [&](const Point& p) {
                return Num::Lt(p.pnt[TI->split.second], TI->split.first);
            }).begin();

        bool putTomb =
            hasTomb &&
            (BT::ImbalanceNode(TI->left->size - (split_iter - In.begin()),
                               TI->size - In.size()) ||
             TI->size - In.size() < BT::kThinLeaveWrap);
        hasTomb = putTomb ? false : hasTomb;
        assert(putTomb ? (!hasTomb) : true);

        // assert(this->_split_rule == MAX_STRETCH_DIM ||
        //        (this->_split_rule == ROTATE_DIM && d == TI->split.second));
        DimsType nextDim = (d + 1) % BT::kDim;

        // TODO: change this to optimized serial version
        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        auto [L, Lbox] = BatchDeleteRecursive(
            TI->left, lbox, In.cut(0, split_iter - In.begin()),
            Out.cut(0, split_iter - In.begin()), nextDim, hasTomb);
        auto [R, Rbox] = BatchDeleteRecursive(
            TI->right, rbox, In.cut(split_iter - In.begin(), n),
            Out.cut(split_iter - In.begin(), n), nextDim, hasTomb);

        TI->aug = hasTomb ? false : TI->size > BT::kSerialBuildCutoff;
        BT::template UpdateInterior<Interior>(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        // NOTE: rebuild
        if (putTomb) {
            assert(TI->size == T->size);
            assert(BT::ImbalanceNode(TI->left->size, TI->size) ||
                   TI->size < BT::kThinLeaveWrap);
            return BT::template RebuildSingleTree<Leaf, Interior, false>(T, d);
        }

        return NodeBox(T, BT::GetBox(Lbox, Rbox));
    }

    typename BT::template InnerTree<Leaf, Interior> IT;
    // IT.init();
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    auto tree_nodes = parlay::sequence<NodeBox>::uninitialized(IT.tags_num);
    auto boxs = parlay::sequence<Box>::uninitialized(IT.tags_num);

    IT.TagInbalanceNodeDeletion(boxs, bx, hasTomb);

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
                IT.tags[IT.rev_tag[i]].first, boxs[i],
                Out.cut(start, start + IT.sums[i]),
                In.cut(start, start + IT.sums[i]), nextDim,
                IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
        },
        1);

    BucketType beatles = 0;
    return DeleteInnerTree(1, IT.tags, tree_nodes, beatles, IT.rev_tag, d);
}

// NOTE: only sieve the Points, without rebuilding the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::BatchDeleteRecursive(
    Node* T, const typename KdTree<Point, SplitRule, kBDO>::Box& bx, Slice In,
    Slice Out, DimsType d, PartialCoverTag) {
    size_t n = In.size();

    if (n == 0) return NodeBox(T, bx);

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);

        if (TL->is_dummy) {  // NOTE: need to check whether all Points are in
                             // the Leaf
            assert(T->is_leaf);

            // PERF: slow when In.size() is large
            for (size_t i = 0; TL->size && i < In.size(); i++) {
                if (TL->pts[0] == In[i]) {
                    TL->size -= 1;
                }
            }
            assert(TL->size >= 0);
            return NodeBox(T, Box(TL->pts[0], TL->pts[0]));
        }

        auto it = TL->pts.begin(), end = TL->pts.begin() + TL->size;
        for (int i = 0; TL->size && i < In.size(); i++) {
            it = std::ranges::find(TL->pts.begin(), end, In[i]);
            if (it != end) {  // NOTE: find a Point
                std::ranges::iter_swap(it, --end);
                TL->size -= 1;
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

        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        auto [L, Lbox] = BatchDeleteRecursive(
            TI->left, lbox, In.cut(0, split_iter - In.begin()),
            Out.cut(0, split_iter - In.begin()), nextDim, PartialCoverTag());
        auto [R, Rbox] = BatchDeleteRecursive(
            TI->right, rbox, In.cut(split_iter - In.begin(), n),
            Out.cut(split_iter - In.begin(), n), nextDim, PartialCoverTag());

        BT::template UpdateInterior<Interior>(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        return NodeBox(T, BT::GetBox(Lbox, Rbox));
    }

    typename BT::template InnerTree<Leaf, Interior> IT;
    // IT.init();
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    auto tree_nodes = parlay::sequence<NodeBox>::uninitialized(IT.tags_num);
    auto boxs = parlay::sequence<Box>::uninitialized(IT.tags_num);

    // NOTE: never set tomb, this equivalent to only calcualte the bounding box,
    IT.TagInbalanceNodeDeletion(boxs, bx, false);

    parlay::parallel_for(
        0, IT.tags_num,
        // NOTE: i is the index of the tags
        [&](size_t i) {
            // assert( IT.sums_tree[IT.rev_tag[i]] == IT.sums[i] );
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
            tree_nodes[i] = BatchDeleteRecursive(
                IT.tags[IT.rev_tag[i]].first, boxs[i],
                Out.cut(start, start + IT.sums[i]),
                In.cut(start, start + IT.sums[i]), nextDim, PartialCoverTag());
        },
        1);

    BucketType beatles = 0;
    return UpdateInnerTree(1, IT.tags, tree_nodes, beatles, IT.rev_tag);
}

}  // namespace cpdd
