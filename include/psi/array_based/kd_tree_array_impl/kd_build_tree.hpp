#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_

#include "array_based/kd_tree_array.h"

namespace psi {
namespace array_based {

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
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), 0, this->tree_box_);
  assert(this->root_ != BT::NULL_INDEX);
  return;
}

template <typename TypeTrait>
auto KdTreeArray<TypeTrait>::BuildRecursive(Slice In, Slice Out, DimsType dim,
                                            Box const& bx) -> NodeIndex {
  // TODO: Implement
  // Recursive build using array indices
  // Similar to pointer-based version but returns NodeIndex

  return BT::NULL_INDEX;
}

template <typename TypeTrait>
auto KdTreeArray<TypeTrait>::SerialBuildRecursive(Slice In, Slice Out,
                                                  DimsType dim, Box const& bx)
    -> NodeIndex {
  // TODO: Implement
  // Serial build for small subtrees

  return BT::NULL_INDEX;
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

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
