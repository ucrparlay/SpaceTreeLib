#ifndef PSTP_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_
#define PSTP_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_

#include "../kd_tree.h"
#include "parlay/slice.h"

namespace pstp {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(Slice A) {
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

// NOTE: return the updated Node
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsertRecursive(
    Node* T, Slice In, Slice Out, DimsType d) {
  size_t n = In.size();

  if (n == 0) return T;

  if (T->is_leaf) {
    Leaf* TL = static_cast<Leaf*>(T);
    if ((!TL->is_dummy && n + T->size <= BT::kLeaveWrap) ||
        (TL->is_dummy &&
         parlay::all_of(In, [&](Point const& p) { return p == TL->pts[0]; }))) {
      return BT::template InsertPoints2Leaf<Leaf>(T, In);
    } else {  // PERF: if a nomarl leaf TL cannot handle more duplicates,
              // leave them here
      return BT::template RebuildWithInsert<Leaf, Interior>(T, In, d);
    }
  }

  if (n <= BT::kSerialBuildCutoff) {
    Interior* TI = static_cast<Interior*>(T);
    std::ranges::subrange _2ndGroup =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        });
    size_t split_pos = static_cast<PointsIter>(_2ndGroup.begin()) -
                       static_cast<PointsIter>(In.begin());

    // NOTE: rebuild
    if (split_rule_.AllowRebuild() &&
        BT::ImbalanceNode(TI->left->size + split_pos, TI->size + n)) {
      return BT::template RebuildWithInsert<Leaf, Interior>(T, In, d);
    }

    // NOTE: continue
    Node *L, *R;
    d = split_rule_.NextDimension(d);
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

  IT.TagInbalanceNode([&]([[maybe_unused]] BucketType idx) -> bool {
    auto const TI = static_cast<Interior*>(IT.tags[idx].first);
    return split_rule_.AllowRebuild() &&
           BT::ImbalanceNode(TI->left->size + IT.sums_tree[idx << 1],
                             TI->size + IT.sums_tree[idx]);
  });
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);

  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t s = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          s += IT.sums_tree[IT.rev_tag[j]];
        }

        DimsType next_dim = d, depth = IT.GetDepthByIndex(IT.rev_tag[i]);
        for (BucketType i = 0; i < depth; i++) {
          next_dim = split_rule_.NextDimension(next_dim);
        }

        if (IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1) {
          // NOTE: continue sieve
          tree_nodes[i] = BatchInsertRecursive(
              IT.tags[IT.rev_tag[i]].first,
              Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
              In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), next_dim);
        } else {  // NOTE: launch rebuild subtree
          assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 2);
          assert(IT.tags[IT.rev_tag[i]].first->size +
                     IT.sums_tree[IT.rev_tag[i]] >=
                 0);

          tree_nodes[i] = BT::template RebuildWithInsert<Leaf, Interior>(
              IT.tags[IT.rev_tag[i]].first,
              Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]), next_dim);
        }
      },
      1);

  return IT.template UpdateInnerTree<InnerTree::kUpdatePointer>(tree_nodes);
}

}  // namespace pstp

#endif  // PSTP_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_
