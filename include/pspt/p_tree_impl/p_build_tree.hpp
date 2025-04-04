#ifndef PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_
#define PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "../p_tree.h"
#include "parlay/utilities.h"
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
  // std::cout << sizeof(A[0]) << "\n" << sizeof(std::ref(A[0])) << "\n";
  // std::cout << A[0] << std::endl;
  // auto o = parlay::sequence<Point>::uninitialized(1);
  // parlay::assign_uninitialized(o[0], std::cref(A[0]));
  // std::cout << o[0] << std::endl;

  // parlay::sequence<CpamPair> entries(n);
  // parlay::parallel_for(0, n, [&](int i) {
  //   // entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
  //   entries[i] = {{space_filling_curve_.Encode(A[i]), i}, std::ref(A[i])};
  //   // entries[i] = {{space_filling_curve_.Encode(A[i]), i}, A[i]};
  //   // entries[i] = {P[i]->id, P[i]};
  // });
  // CpamPair a = {{space_filling_curve_.Encode(A[0]), 0}, A[0]};
  // std::tuple<typename CpamEntry::key_t,
  //            std::reference_wrapper<typename CpamEntry::val_t>>
  //     b = {{space_filling_curve_.Encode(A[0]), 0}, std::ref(A[0])};
  // std::cout << sizeof(a) << " " << sizeof(b) << std::endl;
  // auto a = std::ref(A[0]);
  parlay::internal::timer t("");
  auto entries = parlay::tabulate(n, [&](size_t i) {
    // return {{space_filling_curve_.Encode(A[i]), i}, std::ref(A[i])};
    // return {{0, i}, std::ref(A[i])};
    // return std::make_tuple(std::make_pair(0, i), A[i]);
    // return std::make_tuple(std::make_pair(SplitRule::Encode(A[i]), i), A[i]);
    return std::make_tuple(std::make_pair(SplitRule::Encode(A[i]), i),
                           std::ref(A[i]));
  });
  t.next("tabulate");
  // std::cout << sizeof(entries[0]) << std::endl;
  // static_assert(std::is_same<std::reference_wrapper<Point>,
  //                            typename decltype(entries)::value_type>::value);
  // std::cout << sizeof(entries[0]) << std::endl;
  // zmap m1(entries);
  // auto vals = zmap::values(m1);
  this->cpam_aug_map_ = CpamAugMap(entries);
  t.next("build_cpam_aug_map");
  puts("--------------------");
  // return m1;

  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_BUILD_TREE_HPP_
