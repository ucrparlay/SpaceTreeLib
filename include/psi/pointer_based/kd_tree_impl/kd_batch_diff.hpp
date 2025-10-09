#ifndef PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_
#define PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_

#include "../kd_tree.h"

namespace psi {

template <typename TypeTrait>
template <typename Range>
void KdTree<TypeTrait>::BatchDiff(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDiff_(A);
  return;
}

template <typename TypeTrait>
void KdTree<TypeTrait>::BatchDiff_(Slice A) {
  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  Box box = this->tree_box_;

  DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  std::tie(T, this->tree_box_) =
      BatchDiffRecursive(T, box, A, parlay::make_slice(B), d);

  d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  auto prepare_rebuild_func = [&](Node* T, DimsType d, Box const& box) {
    DimsType new_dim = split_rule_.NextDimension(d);
    BoxCut box_cut(box, static_cast<Interior*>(T)->split, true);
    auto left_args = std::make_pair(new_dim, box_cut.GetFirstBoxCut());
    auto right_args = std::make_pair(new_dim, box_cut.GetSecondBoxCut());
    return std::make_pair(std::move(left_args), std::move(right_args));
  };

  this->root_ = BT::template RebuildTreeRecursive<Leaf, Interior, false>(
      T, prepare_rebuild_func, this->split_rule_.AllowRebuild(), d,
      this->tree_box_);

  return;
}

template <typename TypeTrait>
auto KdTree<TypeTrait>::BatchDiffRecursive(
    Node* T, typename KdTree<TypeTrait>::Box const& box, Slice In, Slice Out,
    DimsType d) -> NodeBox {
  size_t n = In.size();

  if (n == 0) return NodeBox(T, box);

  if (T->is_leaf) {
    return BT::template DiffPoints4Leaf<Leaf, NodeBox>(T, In);
  }

  if (In.size() <= BT::kSerialBuildCutoff) {
    Interior* TI = static_cast<Interior*>(T);
    PointsIter split_iter =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        }).begin();

    DimsType next_dim = split_rule_.NextDimension(d);

    BoxCut box_cut(box, TI->split, true);
    auto [L, Lbox] = BatchDiffRecursive(
        TI->left, box_cut.GetFirstBoxCut(), In.cut(0, split_iter - In.begin()),
        Out.cut(0, split_iter - In.begin()), next_dim);
    auto [R, Rbox] =
        BatchDiffRecursive(TI->right, box_cut.GetSecondBoxCut(),
                           In.cut(split_iter - In.begin(), n),
                           Out.cut(split_iter - In.begin(), n), next_dim);

    bool const force_parallel_flag = TI->size > BT::kSerialBuildCutoff;
    BT::template UpdateInterior<Interior>(T, L, R);
    assert(T->size == L->size + R->size && TI->split.second >= 0 &&
           TI->is_leaf == false);

    if (BT::SparcyNode(0, TI->size) ||
        (split_rule_.AllowRebuild() &&
         BT::ImbalanceNode(TI->left->size, TI->size))) {
      TI->SetParallelFlag(force_parallel_flag);
    }

    return NodeBox(T, Geo::GetBox(Lbox, Rbox));
  }

  InnerTree IT;
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);

  IT.TagOversizedNodes();

  parlay::parallel_for(
      0, IT.tags_num,
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        DimsType next_dim = d, depth = IT.GetDepthByIndex(IT.rev_tag[i]);
        for (BucketType i = 0; i < depth; i++) {
          next_dim = split_rule_.NextDimension(next_dim);
        }

        assert(Geo::WithinBox(
            Geo::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first),
            IT.GetBoxByRegionIdx(IT.rev_tag[i], box)));

        tree_nodes[i] =
            BatchDiffRecursive(IT.tags[IT.rev_tag[i]].first,
                               IT.GetBoxByRegionIdx(IT.rev_tag[i], box),
                               Out.cut(start, start + IT.sums[i]),
                               In.cut(start, start + IT.sums[i]), next_dim);
      },
      1);

  return IT.template UpdateInnerTreePointersWithBox<false>(tree_nodes);
}

}  // namespace psi

#endif  // PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_
