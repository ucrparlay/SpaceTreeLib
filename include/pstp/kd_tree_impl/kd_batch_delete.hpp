#ifndef PSTP_KD_TREE_IMPL_KD_BATCH_DELETE_HPP
#define PSTP_KD_TREE_IMPL_KD_BATCH_DELETE_HPP

#include "../kd_tree.h"

namespace pstp {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDelete_(A);
  return;
}

// NOTE: assume all Points are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete_(Slice A) {
  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  Box bx = this->tree_box_;
  DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  std::tie(this->root_, this->tree_box_) =
      BatchDeleteRecursive(T, bx, A, parlay::make_slice(B), d, 1);
  return;
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
// WARN: the param d can be only used when rotate cutting is applied
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBox
KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDeleteRecursive(
    Node* T,
    typename KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::Box const& bx,
    Slice In, Slice Out, DimsType d, bool has_tomb) {
  size_t n = In.size();

  if (n == 0) {
    assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), bx));
    return NodeBox(T, bx);
  }

  // INFO: it can be used to accelerate the whole deletion process
  if (n == T->size) {
    if (has_tomb) {  // rebuild this subtree
      BT::template DeleteTreeRecursive<Leaf, Interior>(T);
      return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(), BT::GetEmptyBox());
    }
    // within a rebuild tree
    if (!T->is_leaf) {
      auto TI = static_cast<Interior*>(T);
      // WARN: only set the flag for root, the remaining tree is still unset
      TI->SetParallelFlag(T->size > BT::kSerialBuildCutoff);
    }
    T->size = 0;
    return NodeBox(T, bx);
  }

  if (T->is_leaf) {
    return BT::template DeletePoints4Leaf<Leaf, NodeBox>(T, In);
  }

  // if (1) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    Interior* TI = static_cast<Interior*>(T);
    PointsIter split_iter =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        }).begin();

    // NOTE: put the tomb if the remaining points number are below
    // kThinLeaveWrap or inbalance
    bool putTomb =
        BT::SparcyNode(In.size(), TI->size) ||
        (split_rule_.AllowRebuild() && has_tomb &&
         BT::ImbalanceNode(TI->left->size - (split_iter - In.begin()),
                           TI->size - In.size()));
    has_tomb = putTomb ? false : has_tomb;
    assert(putTomb ? (!has_tomb) : true);

    // DimsType next_dim = (d + 1) % BT::kDim;
    DimsType next_dim = split_rule_.NextDimension(d);
    BoxCut box_cut(bx, TI->split, true);

    auto [L, Lbox] = BatchDeleteRecursive(
        TI->left, box_cut.GetFirstBoxCut(), In.cut(0, split_iter - In.begin()),
        Out.cut(0, split_iter - In.begin()), next_dim, has_tomb);
    auto [R, Rbox] = BatchDeleteRecursive(TI->right, box_cut.GetSecondBoxCut(),
                                          In.cut(split_iter - In.begin(), n),
                                          Out.cut(split_iter - In.begin(), n),
                                          next_dim, has_tomb);

    TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
    BT::template UpdateInterior<Interior>(T, L, R);

    assert(T->size == L->size + R->size && TI->split.second >= 0 &&
           TI->is_leaf == false);

    // NOTE: rebuild
    if (putTomb) {
      assert(BT::SparcyNode(0, TI->size) ||
             (split_rule_.AllowRebuild() &&
              BT::ImbalanceNode(TI->left->size, TI->size)));
      auto const new_box = BT::GetBox(Lbox, Rbox);
      assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), new_box));
      return NodeBox(
          BT::template RebuildSingleTree<Leaf, Interior, false>(T, d, new_box),
          new_box);
    }

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
  auto box_seq = parlay::sequence<Box>::uninitialized(IT.tags_num);

  [[maybe_unused]] auto [re_num, tot_re_size] = IT.TagInbalanceNodeDeletion(
      box_seq, bx, has_tomb, [&](BucketType idx) -> bool {
        Interior* TI = static_cast<Interior*>(IT.tags[idx].first);
        return BT::SparcyNode(IT.sums_tree[idx], TI->size) ||
               (split_rule_.AllowRebuild() &&
                BT::ImbalanceNode(TI->left->size - IT.sums_tree[idx << 1],
                                  TI->size - IT.sums_tree[idx]));
      });

  assert(re_num <= IT.tags_num);

  // NOTE: delete the points in the tree
  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t start = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
        assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
        assert(BT::WithinBox(
            BT::GetBox(Out.cut(start, start + IT.sums[i])),
            BT::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first)));

        DimsType next_dim = d, depth = IT.GetDepthByIndex(IT.rev_tag[i]);
        for (BucketType i = 0; i < depth; i++) {
          next_dim = split_rule_.NextDimension(next_dim);
        }

        tree_nodes[i] = BatchDeleteRecursive(
            IT.tags[IT.rev_tag[i]].first, box_seq[i],
            Out.cut(start, start + IT.sums[i]),
            In.cut(start, start + IT.sums[i]), next_dim,
            IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
      },
      1);

  // NOTE: handling of rebuild (in parallel)
  // NOTE: get new box for skeleton root and rebuild nodes
  Box const new_box =
      std::get<1>(IT.template UpdateInnerTree<InnerTree::kTagRebuildNode>(
          tree_nodes, box_seq));
  assert(IT.tags_num == re_num);

  // NOTE: launch the rebuild in parallel
  parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
    assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 3);

    if (IT.tags[IT.rev_tag[i]].first->size == 0) {  // NOTE: empty
      BT::template DeleteTreeRecursive<Leaf, Interior, false>(
          IT.tags[IT.rev_tag[i]].first);
      IT.tags[IT.rev_tag[i]].first = AllocEmptyLeafNode<Slice, Leaf>();
    } else {  // NOTE: rebuild
      assert(BT::WithinBox(
          BT::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first),
          box_seq[i]));

      DimsType next_dim = d, depth = IT.GetDepthByIndex(IT.rev_tag[i]);
      for (BucketType i = 0; i < depth; i++) {
        next_dim = split_rule_.NextDimension(next_dim);
      }
      IT.tags[IT.rev_tag[i]].first =
          BT::template RebuildSingleTree<Leaf, Interior, false>(
              IT.tags[IT.rev_tag[i]].first, next_dim, box_seq[i]);
    }
  });  // PERF: let the parlay decide granularity

  // PARA: op == 0 -> toggle whether under a rebuild tree
  // op == 1 -> query current status
  auto const new_root = std::get<0>(
      IT.template UpdateInnerTree<InnerTree::kPostDelUpdate>(tree_nodes));
  return NodeBox(new_root, new_box);
}

}  // namespace pstp

#endif  // PSTP_KD_TREE_IMPL_KD_BATCH_DELETE_HPP
