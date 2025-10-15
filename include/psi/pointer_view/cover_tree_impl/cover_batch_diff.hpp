#ifndef PSI_POINTER_VIEW_COVER_TREE_IMPL_COVER_BATCH_DIFF_HPP_
#define PSI_POINTER_VIEW_COVER_TREE_IMPL_COVER_BATCH_DIFF_HPP_

#include <tuple>

#include "../cover_tree.h"

namespace psi {
namespace pointer_view {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiff(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  return;
}

}  // namespace pointer_view
}  // namespace psi

#endif  // PSI_POINTER_VIEW_COVER_TREE_IMPL_COVER_BATCH_DIFF_HPP_
