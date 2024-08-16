#pragma once

#include "../kd_tree.h"

namespace cppd {
template<typename Point>
void KdTree<Point>::BatchInsert(Slice A, const DimsType BT::kDim) {
    if (this->root == nullptr) {
        return build(A, BT::kDim);
    }

    Points B = Points::uninitialized(A.size());
    Node* T = this->root;
    /*box b = get_box(A);*/
    this->bbox = get_box(this->bbox, get_box(A));
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    this->root = BatchInsertRecursive(T, A, B.cut(0, A.size()), d, BT::kDim);
    assert(this->root != NULL);
    return;
}

template<typename Point>
typename KdTree<Point>::Node*
KdTree<Point>::UpdateInnerTreeByTag(
    BucketType idx, const NodeTag& tags, parlay::sequence<Node*>& treeNodes,
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

template<typename Point>
typename KdTree<Point>::Node*
KdTree<Point>::RebuildWithInsert(Node* T, Slice In, const DimsType d,
                                           const DimsType BT::kDim) {
    uint_fast8_t curDim = pick_rebuild_dim(T, d, BT::kDim);
    Points wo = Points::uninitialized(T->size + In.size());
    Points wx = Points::uninitialized(T->size + In.size());
    parlay::parallel_for(0, In.size(), [&](size_t j) { wx[j] = In[j]; });
    flatten(T, wx.cut(In.size(), wx.size()));
    delete_tree_recursive(T);
    return build_recursive(parlay::make_slice(wx), parlay::make_slice(wo),
                           curDim, BT::kDim, get_box(parlay::make_slice(wx)));
}

//* return the updated Node
template<typename Point>
typename KdTree<Point>::Node*
KdTree<Point>::BatchInsertRecursive(Node* T, Slice In, Slice Out,
                                            DimsType d, const DimsType BT::kDim) {
    size_t n = In.size();

    if (n == 0) return T;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        if (!TL->is_dummy && n + TL->size <= kLeaveWrap) {
            assert(T->size == TL->size);
            if (TL->pts.size() == 0) {
                TL->pts = Points::uninitialized(kLeaveWrap);
            }
            for (int i = 0; i < n; i++) {
                TL->pts[TL->size + i] = In[i];
            }
            TL->size += n;
            return T;
        } else {
            return RebuildWithInsert(T, In, d, BT::kDim);
        }
    }

    if (n <= BT::SerialBuildCutoff) {
        Interior* TI = static_cast<Interior*>(T);
        auto _2ndGroup = std::ranges::partition(In, [&](const Point& p) {
            return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });

        //* rebuild
        if (inbalance_node(TI->left->size + _2ndGroup.begin() - In.begin(),
                           TI->size + n)) {
            return RebuildWithInsert(T, In, d, BT::kDim);
        }
        //* continue
        Node *L, *R;
        d = (d + 1) % BT::kDim;
        L = BatchInsertRecursive(
            TI->left, In.cut(0, _2ndGroup.begin() - In.begin()),
            Out.cut(0, _2ndGroup.begin() - In.begin()), d, BT::kDim);
        R = BatchInsertRecursive(
            TI->right, In.cut(_2ndGroup.begin() - In.begin(), n),
            Out.cut(_2ndGroup.begin() - In.begin(), n), d, BT::kDim);
        UpdateInterior(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);
        return T;
    }

    //@ assign each Node a tag
    InnerTree IT;
    IT.init();
    assert(IT.rev_tag.size() == BT::kBucketNum);
    IT.assign_node_tag(T, 1);
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);

    seieve_points(In, Out, n, IT.tags, IT.sums, IT.tagsNum);

    IT.tag_inbalance_node();
    assert(IT.tagsNum > 0 && IT.tagsNum <= BT::kBucketNum);
    auto treeNodes = parlay::sequence<Node*>::uninitialized(IT.tagsNum);

    parlay::parallel_for(
        0, IT.tagsNum,
        [&](size_t i) {
            size_t s = 0;
            for (int j = 0; j < i; j++) {
                s += IT.sums_tree[IT.rev_tag[j]];
            }

            DimsType nextDim = (d + IT.get_depth_by_index(IT.rev_tag[i])) % BT::kDim;
            if (IT.tags[IT.rev_tag[i]].second ==
                BT::kBucketNum + 1) {  // NOTE: continue sieve
                treeNodes[i] = BatchInsertRecursive(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
                    In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim, BT::kDim);
            } else {  // NOTE: launch rebuild subtree
                assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 2);
                assert(IT.tags[IT.rev_tag[i]].first->size +
                           IT.sums_tree[IT.rev_tag[i]] >=
                       0);

                treeNodes[i] = RebuildWithInsert(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim, BT::kDim);
            }
        },
        1);

    BucketType beatles = 0;
    return UpdateInnerTreeByTag(1, IT.tags, treeNodes, beatles, IT.rev_tag);
}
// NOTE: update the info of T by new children L and R
template<typename Point>
inline void KdTree<Point>::UpdateInterior(
    typename KdTree<Point>::Node* T,
    typename KdTree<Point>::Node* L,
    typename KdTree<Point>::Node* R) {
    assert(!T->is_leaf);
    Interior* TI = static_cast<Interior*>(T);
    TI->size = L->size + R->size;
    TI->left = L;
    TI->right = R;
    return;
}

// NOTE: retrive the bucket tag of Point p from the skeleton tags
template<typename Point>
uint_fast8_t KdTree<Point>::retrive_tag(const Point& p,
                                                const NodeTag& tags) {
    uint_fast8_t k = 1;
    Interior* TI;
    while (k <= PIVOT_NUM && (!tags[k].first->is_leaf)) {
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
template<typename Point>
void KdTree<Point>::seieve_points(Slice A, Slice B, const size_t n,
                                          const NodeTag& tags,
                                          parlay::sequence<balls_type>& sums,
                                          const BucketType tagsNum) {
    size_t num_block = (n + BLOCK_SIZE - 1) >> LOG2_BASE;
    parlay::sequence<parlay::sequence<balls_type>> offset(
        num_block, parlay::sequence<balls_type>(tagsNum));
    assert(offset.size() == num_block && offset[0].size() == tagsNum &&
           offset[0][0] == 0);
    parlay::parallel_for(0, num_block, [&](size_t i) {
        for (size_t j = i << LOG2_BASE; j < std::min((i + 1) << LOG2_BASE, n);
             j++) {
            offset[i][std::move(retrive_tag(A[j], tags))]++;
        }
    });

    sums = parlay::sequence<balls_type>(tagsNum);
    for (size_t i = 0; i < num_block; i++) {
        auto t = offset[i];
        offset[i] = sums;
        for (int j = 0; j < tagsNum; j++) {
            sums[j] += t[j];
        }
    }

    parlay::parallel_for(0, num_block, [&](size_t i) {
        auto v = parlay::sequence<balls_type>::uninitialized(tagsNum);
        int tot = 0, s_offset = 0;
        for (int k = 0; k < tagsNum - 1; k++) {
            v[k] = tot + offset[i][k];
            tot += sums[k];
            s_offset += offset[i][k];
        }
        v[tagsNum - 1] = tot + ((i << LOG2_BASE) - s_offset);
        for (size_t j = i << LOG2_BASE; j < std::min((i + 1) << LOG2_BASE, n);
             j++) {
            B[v[std::move(retrive_tag(A[j], tags))]++] = A[j];
        }
    });

    return;
}

// NOTE: traverse the skeleton tags recursively and update its children to new
// ones
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::update_inner_tree(BucketType idx, const NodeTag& tags,
                                         parlay::sequence<node_box>& treeNodes,
                                         BucketType& p,
                                         const TagNodes& rev_tag) {
    if (tags[idx].second < BT::kBucketNum) {
        assert(rev_tag[p] == idx);
        return treeNodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != nullptr);
    auto [L, Lbox] = update_inner_tree(idx << 1, tags, treeNodes, p, rev_tag);
    auto [R, Rbox] =
        update_inner_tree(idx << 1 | 1, tags, treeNodes, p, rev_tag);
    UpdateInterior(tags[idx].first, L, R);
    return node_box(tags[idx].first, get_box(Lbox, Rbox));
}

// NOTE: flatten a tree then rebuild upon it
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::rebuild_single_tree(Node* T, const DimsType d,
                                           const DimsType BT::kDim,
                                           const bool granularity) {
    Points wo = Points::uninitialized(T->size);
    Points wx = Points::uninitialized(T->size);
    uint_fast8_t curDim = pick_rebuild_dim(T, d, BT::kDim);
    flatten(T, wx.cut(0, T->size), granularity);
    delete_tree_recursive(T, granularity);
    box bx = get_box(parlay::make_slice(wx));
    Node* o = build_recursive(parlay::make_slice(wx), parlay::make_slice(wo),
                              curDim, BT::kDim, bx);
    return node_box(std::move(o), std::move(bx));
}

// NOTE: traverse the tree in parallel and rebuild the imbalanced subtree
template<typename Point>
typename KdTree<Point>::node_box
KdTree<Point>::rebuild_tree_recursive(Node* T, DimsType d,
                                              const DimsType BT::kDim,
                                              const bool granularity) {
    if (T->is_leaf) {
        return node_box(T, get_box(T));
    }

    Interior* TI = static_cast<Interior*>(T);
    if (inbalance_node(TI->left->size, TI->size)) {
        return rebuild_single_tree(T, d, BT::kDim, granularity);
    }

    Node *L, *R;
    box Lbox, Rbox;
    d = (d + 1) % BT::kDim;
    parlay::par_do_if(
        // NOTE: if granularity is disabled, always traverse the tree in
        // parallel
        (granularity && T->size > BT::SerialBuildCutoff) ||
            (!granularity && TI->aug_flag),
        [&] {
            std::tie(L, Lbox) =
                rebuild_tree_recursive(TI->left, d, BT::kDim, granularity);
        },
        [&] {
            std::tie(R, Rbox) =
                rebuild_tree_recursive(TI->right, d, BT::kDim, granularity);
        });

    UpdateInterior(T, L, R);

    return node_box(T, get_box(Lbox, Rbox));
}

}  // namespace cppd
