#pragma once

#include "../kd_tree.h"
#include "parlay/slice.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchInsert(Slice A) {
    if (this->root_ == nullptr) {  // TODO: may check using explicity tag
        return Build_(A);
    }

    Points B = Points::uninitialized(A.size());
    Node* T = this->root_;
    this->tree_box_ = BT::GetBox(this->tree_box_, BT::GetBox(A));
    DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
    this->root_ = BatchInsertRecursive(T, A, B.cut(0, A.size()), d);
    assert(this->root_ != NULL);
    return;
}

// TODO: remove rev_tag
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::UpdateInnerTreePointer(
    BucketType idx, const NodeTagSeq& tags, parlay::sequence<Node*>& tree_nodes,
    BucketType& p) {
    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
        // assert(rev_tag[p] == idx);
        return tree_nodes[p++];
    }

    // assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != NULL);
    Node *L, *R;
    L = UpdateInnerTreePointer(idx << 1, tags, tree_nodes, p);
    R = UpdateInnerTreePointer(idx << 1 | 1, tags, tree_nodes, p);
    BT::template UpdateInterior<Interior>(tags[idx].first, L, R);
    return tags[idx].first;
}

// template<typename Point, typename SplitRule, uint_fast8_t kBDO>
// Node* KdTree<Point, SplitRule, kBDO>::RebuildWithInsert(Node* T, Slice In,
//                                                         const DimsType d) {
//     DimsType curDim = this->split_rule_.FindRebuildDimension(d);
//     Points wx, wo;
//     BT::template PrepareRebuild<Leaf, Interior>(T, In, wx, wo);
//     return BuildRecursive(parlay::make_slice(wx), parlay::make_slice(wo),
//                           curDim, BT::GetBox(parlay::make_slice(wx)));
// }

// NOTE: return the updated Node
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::BatchInsertRecursive(Node* T, Slice In,
                                                           Slice Out,
                                                           DimsType d) {
    size_t n = In.size();

    if (n == 0) return T;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        if (!TL->is_dummy && n + TL->size <= BT::kLeaveWrap) {
            return BT::template InsertPoints2Leaf<Leaf>(T, In);
        } else {
            return BT::template RebuildWithInsert<Leaf, Interior>(T, In, d);
        }
    }

    if (n <= BT::kSerialBuildCutoff) {
        Interior* TI = static_cast<Interior*>(T);
        std::ranges::subrange _2ndGroup =
            std::ranges::partition(In, [&](const Point& p) {
                return Num::Lt(p.pnt[TI->split.second], TI->split.first);
            });
        size_t split_pos = static_cast<PointsIter>(_2ndGroup.begin()) -
                           static_cast<PointsIter>(In.begin());

        // NOTE: rebuild
        if (BT::ImbalanceNode(TI->left->size + split_pos, TI->size + n)) {
            return BT::template RebuildWithInsert<Leaf, Interior>(T, In, d);
        }

        // NOTE: continue
        Node *L, *R;
        d = (d + 1) % BT::kDim;
        L = BatchInsertRecursive(TI->left, In.cut(0, split_pos),
                                 Out.cut(0, split_pos), d);
        R = BatchInsertRecursive(TI->right, In.cut(split_pos, n),
                                 Out.cut(split_pos, n), d);
        BT::template UpdateInterior<Interior>(T, L, R);
        assert(T->size == L->size + R->size && TI->split.second >= 0 &&
               TI->is_leaf == false);
        return T;
    }

    // NOTE: assign each Node a tag
    InnerTree IT(*this);
    assert(IT.rev_tag.size() == BT::kBucketNum);
    IT.AssignNodeTag(T, 1);
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);

    BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                        IT.tags_num);

    IT.TagInbalanceNode();
    assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);

    parlay::parallel_for(
        0, IT.tags_num,
        [&](size_t i) {
            size_t s = 0;
            for (int j = 0; j < i; j++) {
                s += IT.sums_tree[IT.rev_tag[j]];
            }

            DimsType nextDim =
                (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
            if (IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1) {
                // NOTE: continue sieve
                tree_nodes[i] = BatchInsertRecursive(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
                    In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            } else {  // NOTE: launch rebuild subtree
                assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 2);
                assert(IT.tags[IT.rev_tag[i]].first->size +
                           IT.sums_tree[IT.rev_tag[i]] >=
                       0);

                tree_nodes[i] = BT::template RebuildWithInsert<Leaf, Interior>(
                    IT.tags[IT.rev_tag[i]].first,
                    Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), nextDim);
            }
        },
        1);

    return IT.template UpdateInnerTree<InnerTree::kPointer>(tree_nodes,
                                                            [&]() {});
}

}  // namespace cpdd
