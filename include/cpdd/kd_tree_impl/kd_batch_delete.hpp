#pragma once

#include "../kd_tree.h"

namespace cpdd {
// NOTE: default batch delete
template<typename Point>
void KdTree<Point>::BatchDelete(Slice A, const DimsType BT::kDim) {
    BatchDelete(A, BT::kDim, FullCoveredTag());
    // BatchDelete(A, BT::kDim, PartialCoverTag());
    return;
}

// NOTE: assume all Points are fully covered in the tree
template<typename Point>
void KdTree<Point>::BatchDelete(Slice A, const DimsType BT::kDim,
                                        FullCoveredTag) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root;
    box bx = this->bbox;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    std::tie(this->root, this->bbox) = BatchDeleteRecursive(
        T, bx, A, parlay::make_slice(B), d, BT::kDim, 1, FullCoveredTag());
    return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template<typename Point>
void KdTree<Point>::BatchDelete(Slice A, const DimsType BT::kDim,
                                        PartialCoverTag) {
    Points B = Points::uninitialized(A.size());
    Node* T = this->root;
    box bx = this->bbox;
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    // NOTE: first sieve the Points
    std::tie(T, this->bbox) = BatchDeleteRecursive(
        T, bx, A, parlay::make_slice(B), d, BT::kDim, PartialCoverTag());
    // NOTE: then rebuild the tree with full parallelsim
    std::tie(this->root, bx) = rebuild_tree_recursive(T, d, BT::kDim, false);
    assert(bx == this->bbox);

    return;
}

// NOTE: the Node which needs to be rebuilt has tag BT::kBucketNum+3
// NOTE: the bucket Node whose ancestor has been rebuilt has tag BT::kBucketNum+2
// NOTE: the bucket Node whose ancestor has not been ... has BT::kBucketNum+1
// NOTE: otherwise, it's BT::kBucketNum
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::DeleteInnerTree(BucketType idx, const NodeTag& tags,
                                         parlay::sequence<node_box>& treeNodes,
                                         BucketType& p,
                                         const TagNodes& rev_tag, DimsType d,
                                         const DimsType BT::kDim) {
    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        assert(rev_tag[p] == idx);
        assert(tags[idx].second == BT::kBucketNum + 1 ||
               tags[idx].first->size > BT::SerialBuildCutoff ==
                   static_cast<Interior*>(tags[idx].first)->aug_flag);
        return treeNodes[p++];  // WARN: this blocks the parallelsim
    }

    auto [L, Lbox] = DeleteInnerTree(idx << 1, tags, treeNodes, p, rev_tag,
                                       (d + 1) % BT::kDim, BT::kDim);
    auto [R, Rbox] = DeleteInnerTree(idx << 1 | 1, tags, treeNodes, p,
                                       rev_tag, (d + 1) % BT::kDim, BT::kDim);

    assert(tags[idx].first->size > BT::SerialBuildCutoff ==
           static_cast<Interior*>(tags[idx].first)->aug_flag);
    update_interior(tags[idx].first, L, R);

    if (tags[idx].second == BT::kBucketNum + 3) {  // NOTE: launch rebuild
        Interior const* TI = static_cast<Interior*>(tags[idx].first);
        assert(inbalance_node(TI->left->size, TI->size) ||
               TI->size < THIN_LEAVE_WRAP);

        if (tags[idx].first->size == 0) {  // NOTE: special judge for empty tree
            delete_tree_recursive(tags[idx].first, false);
            return node_box(alloc_empty_leaf(), get_empty_box());
        }

        return rebuild_single_tree(tags[idx].first, d, BT::kDim, false);
    }

    return node_box(tags[idx].first, get_box(Lbox, Rbox));
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree WARN: the param d can be only used when rotate cutting is applied
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::BatchDeleteRecursive(
    Node* T, const typename KdTree<Point>::box& bx, Slice In, Slice Out,
    DimsType d, const DimsType BT::kDim, bool hasTomb, FullCoveredTag) {
    size_t n = In.size();

    if (n == 0) {
        assert(within_box(get_box(T), bx));
        return node_box(T, bx);
    }

    // INFO: may can be used to accelerate the whole deletion process
    // if ( n == T->size ) {
    //     if ( hasTomb ) {
    //         delete_tree_recursive( T );
    //         return node_box( alloc_empty_leaf(), get_empty_box() );
    //     }
    //     T->size = 0;  //* lazy mark
    //     return node_box( T, get_empty_box() );
    // }

    if (T->is_leaf) {
        assert(T->size >= In.size());
        Leaf* TL = static_cast<Leaf*>(T);

        if (TL->is_dummy) {
            assert(T->is_leaf);
            assert(In.size() <=
                   T->size);  // WARN: cannot delete more Points then there are
            T->size -= In.size();  // WARN: this assumes that In\in T
            return node_box(T, box(TL->pts[0], TL->pts[0]));
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
        return node_box(T, get_box(TL->pts.cut(0, TL->size)));
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
             TI->size - In.size() < THIN_LEAVE_WRAP);
        hasTomb = putTomb ? false : hasTomb;
        assert(putTomb ? (!hasTomb) : true);

        assert(this->_split_rule == MAX_STRETCH_DIM ||
               (this->_split_rule == ROTATE_DIM && d == TI->split.second));
        DimsType nextDim = (d + 1) % BT::kDim;

        box lbox(bx), rbox(bx);
        lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
        rbox.first.pnt[TI->split.second] = TI->split.first;

        auto [L, Lbox] = BatchDeleteRecursive(
            TI->left, lbox, In.cut(0, _2ndGroup.begin() - In.begin()),
            Out.cut(0, _2ndGroup.begin() - In.begin()), nextDim, BT::kDim, hasTomb,
            FullCoveredTag());
        auto [R, Rbox] = BatchDeleteRecursive(
            TI->right, rbox, In.cut(_2ndGroup.begin() - In.begin(), n),
            Out.cut(_2ndGroup.begin() - In.begin(), n), nextDim, BT::kDim, hasTomb,
            FullCoveredTag());

        TI->aug_flag = hasTomb ? false : TI->size > this->BT::SerialBuildCutoff;
        update_interior(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);

        // NOTE: rebuild
        if (putTomb) {
            assert(TI->size == T->size);
            assert(inbalance_node(TI->left->size, TI->size) ||
                   TI->size < THIN_LEAVE_WRAP);
            return rebuild_single_tree(T, d, BT::kDim, false);
        }

        return node_box(T, get_box(Lbox, Rbox));
    }

    InnerTree IT;
    IT.init();
    IT.assign_node_tag(T, 1);
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);
    seieve_points(In, Out, n, IT.tags, IT.sums, IT.tagsNum);

    auto treeNodes = parlay::sequence<node_box>::uninitialized(IT.tagsNum);
    auto boxs = parlay::sequence<box>::uninitialized(IT.tagsNum);

    IT.tag_inbalance_node_deletion(boxs, bx, hasTomb);

    parlay::parallel_for(
        0, IT.tagsNum,
        [&](size_t i) {
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
            assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
            assert(within_box(get_box(Out.cut(start, start + IT.sums[i])),
                              get_box(IT.tags[IT.rev_tag[i]].first)));

            DimsType nextDim = (d + IT.get_depth_by_index(IT.rev_tag[i])) % BT::kDim;

            treeNodes[i] = BatchDeleteRecursive(
                IT.tags[IT.rev_tag[i]].first, boxs[i],
                Out.cut(start, start + IT.sums[i]),
                In.cut(start, start + IT.sums[i]), nextDim, BT::kDim,
                IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1,
                FullCoveredTag());
        },
        1);

    BucketType beatles = 0;
    return DeleteInnerTree(1, IT.tags, treeNodes, beatles, IT.rev_tag, d,
                             BT::kDim);
}

// NOTE: only sieve the Points, without rebuilding the tree
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::BatchDeleteRecursive(
    Node* T, const typename KdTree<Point>::box& bx, Slice In, Slice Out,
    DimsType d, const DimsType BT::kDim, PartialCoverTag) {
    size_t n = In.size();

    if (n == 0) return node_box(T, bx);

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
            return node_box(T, box(TL->pts[0], TL->pts[0]));
        }

        auto it = TL->pts.begin(), end = TL->pts.begin() + TL->size;
        for (int i = 0; TL->size && i < In.size(); i++) {
            it = std::ranges::find(TL->pts.begin(), end, In[i]);
            if (it != end) {  // NOTE: find a Point
                std::ranges::iter_swap(it, --end);
                TL->size -= 1;
            }
        }
        return node_box(T, get_box(TL->pts.cut(0, TL->size)));
    }

    if (In.size() <= BT::SerialBuildCutoff) {
        // if (In.size()) {
        Interior* TI = static_cast<Interior*>(T);
        auto _2ndGroup = std::ranges::partition(In, [&](const Point& p) {
            return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });

        DimsType nextDim = (d + 1) % BT::kDim;

        box lbox(bx), rbox(bx);
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

        return node_box(T, get_box(Lbox, Rbox));
    }

    InnerTree IT;
    IT.init();
    IT.assign_node_tag(T, 1);
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);
    seieve_points(In, Out, n, IT.tags, IT.sums, IT.tagsNum);

    auto treeNodes = parlay::sequence<node_box>::uninitialized(IT.tagsNum);
    auto boxs = parlay::sequence<box>::uninitialized(IT.tagsNum);

    // NOTE: never set tomb, this equivalent to only calcualte the bounding box,
    IT.tag_inbalance_node_deletion(boxs, bx, false);

    parlay::parallel_for(
        0, IT.tagsNum,
        // NOTE: i is the index of the tags
        [&](size_t i) {
            // assert( IT.sums_tree[IT.rev_tag[i]] == IT.sums[i] );
            size_t start = 0;
            for (int j = 0; j < i; j++) {
                start += IT.sums[j];
            }

            DimsType nextDim = (d + IT.get_depth_by_index(IT.rev_tag[i])) % BT::kDim;
            treeNodes[i] =
                BatchDeleteRecursive(IT.tags[IT.rev_tag[i]].first, boxs[i],
                                      Out.cut(start, start + IT.sums[i]),
                                      In.cut(start, start + IT.sums[i]),
                                      nextDim, BT::kDim, PartialCoverTag());
        },
        1);

    BucketType beatles = 0;
    return update_inner_tree(1, IT.tags, treeNodes, beatles, IT.rev_tag);
}

}  // namespace cpdd
