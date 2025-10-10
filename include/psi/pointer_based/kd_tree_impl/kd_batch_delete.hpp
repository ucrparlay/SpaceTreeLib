#ifndef PSI_POINTER_BASED_KD_TREE_IMPL_KD_BATCH_DELETE_HPP_
#define PSI_POINTER_BASED_KD_TREE_IMPL_KD_BATCH_DELETE_HPP_

#include "../kd_tree.h"

namespace psi {
namespace pointer_based {

template <typename TypeTrait>
template <typename Range>
void KdTree<TypeTrait>::BatchDelete(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDelete_(A);
  return;
}

template <typename TypeTrait>
void KdTree<TypeTrait>::BatchDelete_(Slice A) {
  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  Box box = this->tree_box_;
  DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  std::tie(this->root_, this->tree_box_) =
      BatchDeleteRecursive(T, box, A, parlay::make_slice(B), d, 1);
  return;
}

template <typename TypeTrait>
auto KdTree<TypeTrait>::BatchDeleteRecursive(
    Node* T, typename KdTree<TypeTrait>::Box const& box, Slice In, Slice Out,
    DimsType d, bool has_tomb) -> NodeBox {
  size_t n = In.size();

  if (n == 0) {
    if constexpr (HasBox<typename Interior::AT> && HasBox<typename Leaf::AT>) {
      assert(Geo::SameBox(Geo::template GetBox<Leaf, Interior>(T),
                          BT::template RetriveBox<Leaf, Interior>(T)));
      return NodeBox(T, BT::template RetriveBox<Leaf, Interior>(T));
    } else {
      assert(Geo::WithinBox(Geo::template GetBox<Leaf, Interior>(T), box));
      return NodeBox(T, box);
    }
  }

  if (n == T->size) {
    if (has_tomb) {
      BT::template DeleteTreeRecursive<Leaf, Interior>(T);
      return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(), Geo::GetEmptyBox());
    }
    if (!T->is_leaf) {
      auto TI = static_cast<Interior*>(T);
      TI->ResetAug();
      TI->SetParallelFlag(T->size > BT::kSerialBuildCutoff);
    } else {
      auto TL = static_cast<Leaf*>(T);
      TL->ResetAug();
    }
    T->size = 0;
    return NodeBox(T, Geo::GetEmptyBox());
  }

  if (T->is_leaf) {
    return BT::template DeletePoints4Leaf<Leaf, NodeBox>(T, In);
  }

  if (In.size() <= BT::kSerialBuildCutoff) {
    Interior* TI = static_cast<Interior*>(T);
    PointsIter split_iter =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        }).begin();

    bool putTomb =
        has_tomb &&
        (BT::SparcyNode(In.size(), TI->size) ||
         (split_rule_.AllowRebuild() &&
          BT::ImbalanceNode(TI->left->size - (split_iter - In.begin()),
                            TI->size - In.size())));
    has_tomb = putTomb ? false : has_tomb;
    assert(putTomb ? (!has_tomb) : true);

    DimsType next_dim = split_rule_.NextDimension(d);
    BoxCut box_cut(box, TI->split, true);

    auto [L, Lbox] = BatchDeleteRecursive(
        TI->left, box_cut.GetFirstBoxCut(), In.cut(0, split_iter - In.begin()),
        Out.cut(0, split_iter - In.begin()), next_dim, has_tomb);
    auto [R, Rbox] = BatchDeleteRecursive(TI->right, box_cut.GetSecondBoxCut(),
                                          In.cut(split_iter - In.begin(), n),
                                          Out.cut(split_iter - In.begin(), n),
                                          next_dim, has_tomb);

    bool const force_parallel_flag = TI->size > BT::kSerialBuildCutoff;

    BT::template UpdateInterior<Interior>(T, L, R);
    TI = static_cast<Interior*>(T);
    assert(T->size == L->size + R->size && TI->split.second >= 0 &&
           TI->is_leaf == false);

    if (!has_tomb) {
      TI->SetParallelFlag(force_parallel_flag);
    }

    if (putTomb) {
      assert(BT::SparcyNode(0, TI->size) ||
             (split_rule_.AllowRebuild() &&
              BT::ImbalanceNode(TI->left->size, TI->size)));
      auto const new_box = Geo::GetBox(Lbox, Rbox);
      assert(Geo::WithinBox(Geo::template GetBox<Leaf, Interior>(T), new_box));
      return NodeBox(this->template RebuildSingleTree<Leaf, Interior, false>(
                         T, d, new_box),
                     new_box);
    }

    return NodeBox(T, Geo::GetBox(Lbox, Rbox));
  }

  InnerTree IT;
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
  auto box_seq = parlay::sequence<Box>::uninitialized(IT.tags_num);

  auto [re_num, tot_re_size] = IT.template TagInbalanceNodeDeletion<true>(
      box_seq, box, has_tomb, [&](BucketType idx) -> bool {
        Interior* TI = static_cast<Interior*>(IT.tags[idx].first);
        return BT::SparcyNode(IT.sums_tree[idx], TI->size) ||
               (split_rule_.AllowRebuild() &&
                BT::ImbalanceNode(TI->left->size - IT.sums_tree[idx << 1],
                                  TI->size - IT.sums_tree[idx]));
      });

  assert(re_num <= IT.tags_num);

  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t start = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
        assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
        assert(Geo::WithinBox(Geo::GetBox(Out.cut(start, start + IT.sums[i])),
                              Geo::template GetBox<Leaf, Interior>(
                                  IT.tags[IT.rev_tag[i]].first)));

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

  Box const new_box = std::get<1>(IT.TagNodesForRebuild(tree_nodes, box_seq));
  assert(IT.tags_num == re_num);

  parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
    assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 3);

    if (IT.tags[IT.rev_tag[i]].first->size == 0) {
      BT::template DeleteTreeRecursive<Leaf, Interior, false>(
          IT.tags[IT.rev_tag[i]].first);
      IT.tags[IT.rev_tag[i]].first = AllocEmptyLeafNode<Slice, Leaf>();
    } else {
      assert(Geo::WithinBox(
          Geo::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first),
          box_seq[i]));

      DimsType next_dim = d, depth = IT.GetDepthByIndex(IT.rev_tag[i]);
      for (BucketType i = 0; i < depth; i++) {
        next_dim = split_rule_.NextDimension(next_dim);
      }
      IT.tags[IT.rev_tag[i]].first =
          this->template RebuildSingleTree<Leaf, Interior, false>(
              IT.tags[IT.rev_tag[i]].first, next_dim, box_seq[i]);
    }
  });

  auto const new_root = std::get<0>(IT.UpdateAfterDeletion(tree_nodes));
  return NodeBox(new_root, new_box);
}

}  // namespace pointer_based
}  // namespace psi

#endif  // PSI_POINTER_BASED_KD_TREE_IMPL_KD_BATCH_DELETE_HPP_