#ifndef PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_
#define PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "../p_tree.h"
#include "parlay/utilities.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::build(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { build_(A); });
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::build_(slice_type A)
{
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
	//   entries[i] = {{space_filling_curve_.encode(A[i]), i},
	//   std::ref(A[i])};
	//   // entries[i] = {{space_filling_curve_.encode(A[i]), i}, A[i]};
	//   // entries[i] = {P[i]->id, P[i]};
	// });
	// CpamPair a = {{space_filling_curve_.encode(A[0]), 0}, A[0]};
	// std::tuple<typename cpam_entry::key_t,
	//            std::reference_wrapper<typename cpam_entry::val_t>>
	//     b = {{space_filling_curve_.encode(A[0]), 0}, std::ref(A[0])};
	// std::cout << sizeof(a) << " " << sizeof(b) << std::endl;
	// auto a = std::ref(A[0]);

	// parlay::internal::timer t("");
	// auto entries = parlay::tabulate(n, [&](size_t i) {
	//   // return {{space_filling_curve_.encode(A[i]), i}, std::ref(A[i])};
	//   // return {{0, i}, std::ref(A[i])};
	//   // return std::make_tuple(std::make_pair(0, i), A[i]);
	//   // return std::make_tuple(std::make_pair(SplitRule::encode(A[i]),
	//   i), A[i]); return std::make_tuple(
	//       std::make_pair(SplitRule::encode(A[i]),
	//       static_cast<id_type>(i)), std::ref(A[i]));
	// });
	// t.next("make_entries");

	// std::cout << sizeof(entries[0]) << std::endl;
	// static_assert(std::is_same<std::reference_wrapper<Point>,
	//                            typename
	//                            decltype(entries)::value_type>::value);
	// std::cout << sizeof(entries[0]) << std::endl;
	// zmap m1(entries);
	// auto vals = zmap::values(m1);
	// parlay::parallel_for(
	//     0, n, [&](size_t i) { A[i].get_aug().code =
	//     SplitRule::encode(A[i]); });
	// this->cpam_aug_map_ = cpam_aug_map_type(A);
	this->cpam_aug_map_ = std::move(cpam_aug_map_type(A));
	// return m1;

	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_
