#ifndef PSTP_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_
#define PSTP_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_

#include "../orth_tree.h"

namespace pstp {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
template <typename Range>
void OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::BatchDelete(
    Range&& In) {
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
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
void OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::BatchDelete_(
    Slice A) {
  Points B = Points::uninitialized(A.size());
  this->root_ = BatchDeleteRecursive(this->root_, A, parlay::make_slice(B),
                                     this->tree_box_, 1);
  return;
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
Node* OrthTree<Point, SplitRule, kMD, kSkHeight,
               kImbaRatio>::BatchDeleteRecursive(Node* T, Slice In, Slice Out,
                                                 Box const& box,
                                                 bool has_tomb) {
  size_t n = In.size();

  if (n == 0) {
    return T;
  }

  // INFO: delete the whole tree directly if its size equals to the input size,
  // can be used to accelerate the whole deletion process
  if (n == T->size) {
    if (has_tomb) {
      BT::template DeleteTreeRecursive<Leaf, Interior>(T);
      return AllocEmptyLeafNode<Slice, Leaf>();
    }
    if (!T->is_leaf) {
      auto TI = static_cast<Interior*>(T);
      TI->SetParallelFlag(T->size > BT::kSerialBuildCutoff);
    }
    T->size = 0;
    return T;
  }

  if (T->is_leaf) {
    return BT::template DeletePoints4Leaf<Leaf, Node*>(T, In);
  }

  if (In.size() <= BT::kSerialBuildCutoff) {
    parlay::sequence<BallsType> sums(kNodeRegions, 0);
    SerialSplitSkeleton(T, In, 0, 1, sums);
    assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

    bool putTomb = has_tomb && (BT::SparcyNode(In.size(), T->size));
    has_tomb = putTomb ? false : has_tomb;
    assert(putTomb ? (!has_tomb) : true);

    auto TI = static_cast<Interior*>(T);
    OrthNodeArr new_nodes;

    size_t start = 0;
    for (DimsType i = 0; i < kNodeRegions; ++i) {
      new_nodes[i] = BatchDeleteRecursive(
          TI->tree_nodes[i], In.cut(start, start + sums[i]),
          Out.cut(start, start + sums[i]), TI->GetBoxByRegionId(i, box),
          has_tomb);
      start += sums[i];
    }

    bool const force_parallel_flag = TI->size > BT::kSerialBuildCutoff;

    BT::template UpdateInterior<Interior>(T, new_nodes);
    assert(T->is_leaf == false);

    if (!has_tomb) {
      TI->SetParallelFlag(force_parallel_flag);
    }

    if (putTomb) {
      // NOTE: the box is the one that passed from the top, which should be the
      // correct one associated with this node
      assert(T->size <= BT::kLeaveWrap);
      assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), box));
      return BT::template RebuildSingleTree<Leaf, Interior, false>(T, box);
    }
    return T;
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  BoxSeq box_seq(IT.tags_num);  // PARA: the box for bucket nodes
  [[maybe_unused]] auto [re_num, tot_re_size] =
      IT.template TagInbalanceNodeDeletion<true>(
          box_seq, box, has_tomb, [&](BucketType idx) -> bool {
            // NOTE: only the sparcy node will be rebuilt
            return BT::SparcyNode(IT.sums_tree[idx], IT.tags[idx].first->size);
          });

  assert(re_num <= IT.tags_num);

  // NOTE: continue seieving points to leaf first
  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);
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

        tree_nodes[i] = BatchDeleteRecursive(
            IT.tags[IT.rev_tag[i]].first, Out.cut(start, start + IT.sums[i]),
            In.cut(start, start + IT.sums[i]), box_seq[i],
            IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
      },
      1);

  // NOTE: handling of rebuild
  // WARN: the rebuild node is on top
  // NOTE: retag the inba-nodes and save the bounding boxes
  [[maybe_unused]] Node* new_node =
      IT.template UpdateInnerTree<InnerTree::kTagRebuildNode>(tree_nodes);
  assert(IT.tags_num == re_num);

  parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
    assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 3);

    if (IT.tags[IT.rev_tag[i]].first->size == 0) {  // NOTE: empty
      BT::template DeleteTreeRecursive<Leaf, Interior, false>(
          IT.tags[IT.rev_tag[i]].first);
      IT.tags[IT.rev_tag[i]].first = AllocEmptyLeafNode<Slice, Leaf>();
    } else {  // NOTE: rebuild
      assert(BT::WithinBox(
          BT::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first),
          IT.GetBoxByRegionIdx(IT.rev_tag[i], box)));
      IT.tags[IT.rev_tag[i]].first =
          BT::template RebuildSingleTree<Leaf, Interior, false>(
              IT.tags[IT.rev_tag[i]].first,
              IT.GetBoxByRegionIdx(IT.rev_tag[i], box));
    }
  });  // PERF: allow the parlay decide the granularity to accelerate the
       // small tree rebuild

  return IT.template UpdateInnerTree<InnerTree::kPostDelUpdate>(tree_nodes);
}

}  // namespace pstp

#endif  // PSTP_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_
