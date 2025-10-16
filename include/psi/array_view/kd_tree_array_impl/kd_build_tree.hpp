#ifndef PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
#define PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_

#include "array_view/kd_tree_array.h"

namespace psi {
namespace array_view {

//==============================================================================
// BUILD OPERATIONS
//==============================================================================

template <typename TypeTrait>
template <typename Range>
void KdTreeArray<TypeTrait>::Build(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A);
}

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = Geo::GetBox(A);
  BuildRecursive(A, B.cut(0, A.size()), 0, this->tree_box_);
  assert(this->root_ != BT::NULL_INDEX);
  return;
}

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::BuildRecursive(Slice In, Slice Out, DimsType dim,
                                            Box const& bx) {
  // TODO: Implement
  // Recursive build using array indices
  // Similar to pointer-based version but returns NodeIndex
}

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::SerialBuildRecursive(Slice In, Slice Out,
                                                  DimsType dim, Box const& box,
                                                  NodeIndex inner_offset,
                                                  NodeIndex leaf_offset) {
  size_t n = In.size();

  if (n == 0) {
    return AnnoteEmptyLeaf(inner_tree_seq_, inner_offset, leaf_offset);
  }

  if (n <= BT::kLeaveWrap) {
    return AnnoteLeaf(inner_tree_seq_, leaf_seq_, In, inner_offset, leaf_offset,
                      BT::kLeaveWrap);
  }

  DimsType next_dim = split_rule_.FindCuttingDimension(box, dim);
  auto [split_iter, split] = split_rule_.SplitInput(In, next_dim, box);

  if (!split.has_value()) {
    if (In.end() == std::ranges::find_if_not(In, [&](Point const& p) {
          return p.SameDimension(In[0]);
        })) {  // dimensions for all points are identical
      if constexpr (IsAugPoint<Point> && Point::IsNonTrivialAugmentation()) {
        return AnnoteLeaf(
            inner_tree_seq_, leaf_seq_, In, inner_offset, leaf_offset,
            std::max(In.size(), static_cast<size_t>(BT::kLeaveWrap)));
      } else {  // the points is identified by only the coordinates
        return AnnoteDummyLeaf(inner_tree_seq_, leaf_seq_, In, inner_offset,
                               leaf_offset);
      }
    } else {  // current dimension is same for all points
      return split_rule_.HandlingUndivide(*this, In, Out, box, dim);
    }
  }

  assert(std::ranges::all_of(In.begin(), split_iter, [&](Point& p) {
    return Num::Lt(p.pnt[split.value().second], split.value().first);
  }));
  assert(std::ranges::all_of(split_iter, In.end(), [&](Point& p) {
    return Num::Geq(p.pnt[split.value().second], split.value().first);
  }));

  BoxCut box_cut(box, split.value(), true);

  next_dim = split_rule_.NextDimension(next_dim);

  auto split_pos = split_iter - In.begin();
  SerialBuildRecursive(In.cut(0, split_pos), Out.cut(0, split_pos), next_dim,
                       box_cut.GetFirstBoxCut(), inner_offset * 2, leaf_offset);
  SerialBuildRecursive(In.cut(split_pos, n), Out.cut(split_pos, n), next_dim,
                       box_cut.GetSecondBoxCut(), inner_offset * 2 + 1,
                       leaf_offset + split_pos);
  AnnoteInterior(inner_tree_seq_, inner_offset, leaf_offset, In.size(),
                 split.value());
  return;
}

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::PickPivots(Slice In, size_t const& n,
                                        SplitterSeq& pivots, DimsType const dim,
                                        BoxSeq& box_seq, Box const& bx) {
  // TODO: Implement
  // Pivot selection logic (can reuse from pointer-based)
}

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::DivideRotate(Slice In, SplitterSeq& pivots,
                                          DimsType dim, BucketType idx,
                                          BoxSeq& box_seq, Box const& bx) {
  // TODO: Implement
  // Partitioning logic (can reuse from pointer-based)
}

}  // namespace array_view
}  // namespace psi

#endif  // PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
