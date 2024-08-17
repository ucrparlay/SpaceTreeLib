#pragma once

#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

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
// NOTE: default batch delete
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Slice A) {
    BatchDelete(A, BT::kDim, FullCoveredTag());
    // BatchDelete(A, BT::kDim, PartialCoverTag());
    return;
}

// NOTE: assume all Points are fully covered in the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Slice A, FullCoveredTag) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    Box bx = this->tree_box_;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    std::tie(this->root_, this->tree_box_) = BatchDeleteRecursive(
        T, bx, A, parlay::make_slice(B), d, 1, FullCoveredTag());
    return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Slice A, PartialCoverTag) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root;
    Box bx = this->bbox;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    // NOTE: first sieve the Points
    std::tie(T, this->bbox) = BatchDeleteRecursive(
        T, bx, A, parlay::make_slice(B), d, PartialCoverTag());
    // NOTE: then rebuild the tree with full parallelsim
    std::tie(this->root, bx) = RebuildTreeRecursive(T, d, false);
    assert(bx == this->bbox);

    return;
}

// NOTE: the Node which needs to be rebuilt has tag BT::kBucketNum+3
// NOTE: the bucket Node whose ancestor has been rebuilt has tag
// BT::kBucketNum+2 NOTE: the bucket Node whose ancestor has not been ... has
// BT::kBucketNum+1 NOTE: otherwise, it's BT::kBucketNum
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
               tags[idx].first->size > BT::SerialBuildCutoff ==
                   static_cast<Interior*>(tags[idx].first)->aug_flag);
        return tree_nodes[p++];  // WARN: this blocks the parallelsim
    }

    auto [L, Lbox] = DeleteInnerTree(idx << 1, tags, tree_nodes, p, rev_tag,
                                     (d + 1) % BT::kDim, BT::kDim);
    auto [R, Rbox] = DeleteInnerTree(idx << 1 | 1, tags, tree_nodes, p, rev_tag,
                                     (d + 1) % BT::kDim, BT::kDim);

    assert(tags[idx].first->size > BT::SerialBuildCutoff ==
           static_cast<Interior*>(tags[idx].first)->aug_flag);
    update_interior(tags[idx].first, L, R);

    if (tags[idx].second == BT::kBucketNum + 3) {  // NOTE: launch rebuild
        Interior const* TI = static_cast<Interior*>(tags[idx].first);
        assert(inbalance_node(TI->left->size, TI->size) ||
               TI->size < BT::kThinLeaveWrap);

        if (tags[idx].first->size == 0) {  // NOTE: special judge for empty tree
            delete_tree_recursive(tags[idx].first, false);
            return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(),
                           BT::GetEmptyBox());
        }

        return rebuild_single_tree(tags[idx].first, d, BT::kDim, false);
    }

    return NodeBox(tags[idx].first, BT::GetBox(Lbox, Rbox));
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree WARN: the param d can be only used when rotate cutting is applied
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::BatchDeleteRecursive(
    Node* T, const typename KdTree<Point, SplitRule, kBDO>::Box& bx, Slice In,
    Slice Out, DimsType d, bool hasTomb, FullCoveredTag) {
    size_t n = In.size();

    if (n == 0) {
        assert(within_box(BT::GetBox(T), bx));
        return NodeBox(T, bx);
    }

    // INFO: may can be used to accelerate the whole deletion process
    // if ( n == T->size ) {
    //     if ( hasTomb ) {
    //         delete_tree_recursive( T );
    //         return NodeBox( alloc_empty_leaf(), get_empty_box() );
    //     }
    //     T->size = 0;  //* lazy mark
    //     return NodeBox( T, get_empty_box() );
    // }

    if (T->is_leaf) {
        assert(T->size >= In.size());
        Leaf* TL = static_cast<Leaf*>(T);

        if (TL->is_dummy) {
            assert(T->is_leaf);
            assert(In.size() <=
                   T->size);  // WARN: cannot delete more Points then there are
            T->size -= In.size();  // WARN: this assumes that In\in T
            return NodeBox(T, box(TL->pts[0], TL->pts[0]));
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

    if (In.size() <= BT::SerialBuildCutoff) {
        Interior* TI = static_cast<Interior*>(T);
        auto _2ndGroup = std::ranges::partition(In, [&](const Point& p) {
            return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });

        bool putTomb =
            hasTomb &&
            (inbalance_node(TI->left->size - (_2ndGroup.begin() - In.begin()),
                            TI->size - In.size()) ||
             TI->size - In.size() < BT::kThinLeaveWrap);
        hasTomb = putTomb ? false : hasTomb;
        assert(putTomb ? (!hasTomb) : true);

        // assert(this->_split_rule == MAX_STRETCH_DIM ||
        //        (this->_split_rule == ROTATE_DIM && d == TI->split.second));
        DimsType nextDim = (d + 1) % BT::kDim;

        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        auto [L, Lbox] = BatchDeleteRecursive(
            TI->left, lbox, In.cut(0, _2ndGroup.begin() - In.begin()),
            Out.cut(0, _2ndGroup.begin() - In.begin()), nextDim, hasTomb,
            FullCoveredTag());
        auto [R, Rbox] = BatchDeleteRecursive(
            TI->right, rbox, In.cut(_2ndGroup.begin() - In.begin(), n),
            Out.cut(_2ndGroup.begin() - In.begin(), n), nextDim, hasTomb,
            FullCoveredTag());

        TI->aug_flag = hasTomb ? false : TI->size > this->BT::SerialBuildCutoff;
        update_interior(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        // NOTE: rebuild
        if (putTomb) {
            assert(TI->size == T->size);
            assert(inbalance_node(TI->left->size, TI->size) ||
                   TI->size < BT::kThinLeaveWrap);
            return rebuild_single_tree(T, d, BT::kDim, false);
        }

        return NodeBox(T, BT::GetBox(Lbox, Rbox));
    }

    typename BT::template InnerTree<Leaf, Interior> IT;
    IT.init();
    IT.assign_node_tag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    seieve_points(In, Out, n, IT.tags, IT.sums, IT.tags_num);

    auto tree_nodes = parlay::sequence<NodeBox>::uninitialized(IT.tags_num);
    auto boxs = parlay::sequence<Box>::uninitialized(IT.tags_num);

    IT.tag_inbalance_node_deletion(boxs, bx, hasTomb);

    parlay::parallel_for(
        0, IT.tags_num,
        [&](size_t i) {
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
            assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
            assert(within_box(BT::GetBox(Out.cut(start, start + IT.sums[i])),
                              BT::GetBox(IT.tags[IT.rev_tag[i]].first)));

            DimsType nextDim =
                (d + IT.get_depth_by_index(IT.rev_tag[i])) % BT::kDim;

            tree_nodes[i] = BatchDeleteRecursive(
                IT.tags[IT.rev_tag[i]].first, boxs[i],
                Out.cut(start, start + IT.sums[i]),
                In.cut(start, start + IT.sums[i]), nextDim, BT::kDim,
                IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1,
                FullCoveredTag());
        },
        1);

    BucketType beatles = 0;
    return DeleteInnerTree(1, IT.tags, tree_nodes, beatles, IT.rev_tag, d,
                           BT::kDim);
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
            return NodeBox(T, box(TL->pts[0], TL->pts[0]));
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

    if (In.size() <= BT::SerialBuildCutoff) {
        // if (In.size()) {
        Interior* TI = static_cast<Interior*>(T);
        auto _2ndGroup = std::ranges::partition(In, [&](const Point& p) {
            return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });

        DimsType nextDim = (d + 1) % BT::kDim;

        Box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        auto [L, Lbox] = BatchDeleteRecursive(
            TI->left, lbox, In.cut(0, _2ndGroup.begin() - In.begin()),
            Out.cut(0, _2ndGroup.begin() - In.begin()), nextDim, BT::kDim,
            PartialCoverTag());
        auto [R, Rbox] = BatchDeleteRecursive(
            TI->right, rbox, In.cut(_2ndGroup.begin() - In.begin(), n),
            Out.cut(_2ndGroup.begin() - In.begin(), n), nextDim, BT::kDim,
            PartialCoverTag());

        update_interior(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        return NodeBox(T, BT::GetBox(Lbox, Rbox));
    }

    typename BT::template InnerTree<Leaf, Interior> IT;
    IT.init();
    IT.assign_node_tag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    seieve_points(In, Out, n, IT.tags, IT.sums, IT.tags_num);

    auto tree_nodes = parlay::sequence<NodeBox>::uninitialized(IT.tags_num);
    auto boxs = parlay::sequence<Box>::uninitialized(IT.tags_num);

    // NOTE: never set tomb, this equivalent to only calcualte the bounding box,
    IT.tag_inbalance_node_deletion(boxs, bx, false);

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
                (d + IT.get_depth_by_index(IT.rev_tag[i])) % BT::kDim;
            tree_nodes[i] =
                BatchDeleteRecursive(IT.tags[IT.rev_tag[i]].first, boxs[i],
                                     Out.cut(start, start + IT.sums[i]),
                                     In.cut(start, start + IT.sums[i]), nextDim,
                                     BT::kDim, PartialCoverTag());
        },
        1);

    BucketType beatles = 0;
    return update_inner_tree(1, IT.tags, tree_nodes, beatles, IT.rev_tag);
}

}  // namespace cpdd
