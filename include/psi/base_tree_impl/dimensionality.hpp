#ifndef PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
#define PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_

#include "../base_tree.h"

namespace psi
{

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline size_t
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_imbalance_ratio()
{
	return static_cast<size_t>(imbalance_ratio);
	// if (auto const env_p = std::getenv("imbalance_ratio")) {
	//   return static_cast<size_t>(std::stoi(env_p));
	// } else {
	//   return static_cast<size_t>(imbalance_ratio);
	// }
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::imbalance_node(
	size_t const l, size_t const n)
{
	if (n == 0)
		return true;
	/*
	 * Compare in double. The old form cast to size_t first, so an
	 * imbalance in (30, 31) truncated to 30 and never fired, and it went
	 * through num_comparator<coord_type>, which exists for coordinates and
	 * carries an epsilon when coord_type is floating point.
	 */
	double const off =
		100.0 * static_cast<double>(l) / static_cast<double>(n) - 50.0;
	return std::abs(off) > static_cast<double>(get_imbalance_ratio());
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::sparse_node(size_t const rm,
								size_t const n)
{
	// PERF: to avoid the case that the new leaf is about 32 and then next
	// slight larger insert will break the leaf
	return n - rm < sparse_leaf_threshold;
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
