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
                                                  DimsType dim,
                                                  Box const& box) {
  size_t n = In.size();

  if (n == 0) return AllocEmptyLeafNode<Slice, Leaf>();

  if (n <= BT::kLeaveWrap) return AllocNormalLeafNode<Slice, Leaf>(In);

  DimsType d = split_rule_.FindCuttingDimension(box, dim);
  auto [split_iter, split] = split_rule_.SplitInput(In, d, box);

  if (!split.has_value()) {
    if (In.end() == std::ranges::find_if_not(In, [&](Point const& p) {
          return p.SameDimension(In[0]);
        })) {
      if constexpr (IsAugPoint<Point>) {
        if constexpr (Point::IsNonTrivialAugmentation()) {
          return AllocFixSizeLeafNode<Slice, Leaf>(
              In, std::max(In.size(), static_cast<size_t>(BT::kLeaveWrap)));
        } else {
          return AllocDummyLeafNode<Slice, Leaf>(In);
        }
      } else {
        return AllocDummyLeafNode<Slice, Leaf>(In);
      }
    } else {
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

  d = split_rule_.NextDimension(d);

  SerialBuildRecursive(In.cut(0, split_iter - In.begin()),
                       Out.cut(0, split_iter - In.begin()), d,
                       box_cut.GetFirstBoxCut());
  SerialBuildRecursive(In.cut(split_iter - In.begin(), n),
                       Out.cut(split_iter - In.begin(), n), d,
                       box_cut.GetSecondBoxCut());
  // return AllocInteriorNode<Interior>(L, R, split.value());
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
