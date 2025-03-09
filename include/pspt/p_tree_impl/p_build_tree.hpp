#ifndef PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_
#define PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "../p_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(Slice A) {
  size_t n = A.size();

  // parlay::parallel_for(0, n, [&](size_t i) {
  //   P[i].morton_id = uRse_hilbert ? P[i].overlap_bits() :
  //   P[i].interleave_bits();
  // });

  parlay::sequence<par> entries(n);
  parlay::parallel_for(0, n, [&](int i) {
    // entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
    entries[i] = {{space_filling_curve_.Encode(A[i]), i}, A[i]};
    // entries[i] = {P[i]->id, P[i]};
  });
  // zmap m1(entries);
  // auto vals = zmap::values(m1);
  this->cpam_aug_map_ = CpamAugMap(entries);
  // return m1;

  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_
