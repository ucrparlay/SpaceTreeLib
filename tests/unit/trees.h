#ifndef PSI_TESTS_UNIT_TREES_H_
#define PSI_TESTS_UNIT_TREES_H_

/*
 * The three trees under one interface, so a property is written once and run
 * against all of them. p_tree needs an augmented point carrying a curve code,
 * which is why the point type is per-kind rather than shared.
 */

#include <parlay/sequence.h>

#include <ostream>
#include <vector>

#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"
#include "psi/orth_tree.h"
#include "psi/p_tree.h"
#include "tests/unit/reference.h"

namespace psi_test
{

template <typename Coord, uint_fast8_t Dim>
struct kd_kind {
	using point_type = psi::basic_point<Coord, Dim>;
	using split_rule =
		psi::orthogonal_split_rule<psi::rotate_dim<point_type>,
					   psi::spatial_median<point_type>>;
	using tree_type =
		psi::kd_tree<psi::tree_traits<point_type, split_rule>>;
	static constexpr char const *name = "kd_tree";
};

template <typename Coord, uint_fast8_t Dim>
struct orth_kind {
	using point_type = psi::basic_point<Coord, Dim>;
	using split_rule =
		psi::orthogonal_split_rule<psi::rotate_dim<point_type>,
					   psi::spatial_median<point_type>>;
	using tree_type =
		psi::orth_tree<psi::orth_tree_traits<point_type, split_rule>>;
	static constexpr char const *name = "orth_tree";
};

/* p_tree keys points by their curve code; the id breaks ties and is what
 * point equality means for it. */
struct curve_code {
	using id_type = int_fast32_t;
	using curve_code_type = uint64_t;

	curve_code() : code(0), id(0)
	{
	}

	void set_member(curve_code_type const &val)
	{
		code = val;
	}

	bool operator<(curve_code const &rhs) const
	{
		return code == rhs.code ? id < rhs.id : code < rhs.code;
	}

	bool operator==(curve_code const &rhs) const
	{
		return id == rhs.id;
	}

	friend std::ostream &operator<<(std::ostream &os, curve_code const &rhs)
	{
		return os << rhs.code << " " << rhs.id;
	}

	curve_code_type code;
	id_type id;
};

template <typename Coord, uint_fast8_t Dim>
struct p_kind {
	using point_type = psi::aug_point<Coord, Dim, curve_code>;
	using split_rule =
		psi::spatial_filling_curve<psi::morton_curve<point_type>>;
	using tree_type = psi::p_tree<psi::tree_traits<point_type, split_rule>>;
	static constexpr char const *name = "p_tree";
};

/*
 * The rules the three kinds above never reach. Both pairs already ship in
 * example/, so these are supported configurations rather than new ones:
 * example/kd_tree.h builds with max_stretch_dim and object_median, and
 * example/p_tree.h with the hilbert curve.
 */
template <typename Coord, uint_fast8_t Dim>
struct kd_stretch_kind {
	using point_type = psi::basic_point<Coord, Dim>;
	using split_rule =
		psi::orthogonal_split_rule<psi::max_stretch_dim<point_type>,
					   psi::object_median<point_type>>;
	using tree_type =
		psi::kd_tree<psi::tree_traits<point_type, split_rule>>;
	static constexpr char const *name = "kd_tree/max_stretch";
};

template <typename Coord, uint_fast8_t Dim>
struct p_hilbert_kind {
	using point_type = psi::aug_point<Coord, Dim, curve_code>;
	using split_rule =
		psi::spatial_filling_curve<psi::hilbert_curve<point_type>>;
	using tree_type = psi::p_tree<psi::tree_traits<point_type, split_rule>>;
	static constexpr char const *name = "p_tree/hilbert";
};

/* Distinct ids matter for p_tree and are harmless elsewhere. */
template <typename Kind>
parlay::sequence<typename Kind::point_type>
make_points(std::vector<typename Kind::point_type> const &src)
{
	parlay::sequence<typename Kind::point_type> out(src.size());
	for (size_t i = 0; i < src.size(); i++) {
		out[i] = src[i];
		if constexpr (requires { out[i].aug.id; }) {
			out[i].aug.id = static_cast<decltype(out[i].aug.id)>(i);
		}
	}
	return out;
}

template <typename Kind>
std::vector<typename Kind::point_type>
gen(size_t n, typename Kind::point_type::coord_type lo,
    typename Kind::point_type::coord_type hi, uint64_t seed)
{
	auto v = random_points<typename Kind::point_type>(n, lo, hi, seed);
	for (size_t i = 0; i < v.size(); i++) {
		if constexpr (requires { v[i].aug.id; }) {
			v[i].aug.id = static_cast<decltype(v[i].aug.id)>(i);
		}
	}
	return v;
}

} // namespace psi_test

#endif /* PSI_TESTS_UNIT_TREES_H_ */
