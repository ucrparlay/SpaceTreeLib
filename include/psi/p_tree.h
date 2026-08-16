#ifndef PSI_P_TREE_H
#define PSI_P_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "base_tree.h"

/*
 * cpam is third party. Silence its noise, but not -Wreturn-local-addr or
 * -Wuninitialized: both currently fire in augmented_node.h and both are real.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "dependence/cpam/cpam.h"
#pragma GCC diagnostic pop

namespace psi
{

template <typename Point, typename SplitRule, uint_fast8_t SkHeight = 6,
	  uint_fast8_t ImbaRatio = 30>
class p_tree
    : public base_tree<Point, p_tree<Point, SplitRule, SkHeight, ImbaRatio>,
		       SkHeight, ImbaRatio>
{
public:
	using base_type =
		base_tree<Point, p_tree<Point, SplitRule, SkHeight, ImbaRatio>,
			  SkHeight, ImbaRatio>;

	using bucket_type = typename base_type::bucket_type;
	using balls_type = typename base_type::balls_type;
	using dims_type = typename base_type::dims_type;
	using bucket_seq_type = typename base_type::bucket_seq_type;
	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	using num_type = num_comparator<coord_type>;
	using slice_type = typename base_type::slice_type;
	using points_type = typename base_type::points_type;
	using points_iter_type = typename base_type::points_iter_type;
	using box_type = typename base_type::box_type;
	using box_seq_type = typename base_type::box_seq_type;

	using hyper_plane_type = typename base_type::hyper_plane_type;
	using hyper_plane_seq_type = typename base_type::hyper_plane_seq_type;
	using node_box_type = typename base_type::node_box_type;
	using node_box_seq_type = typename base_type::node_box_seq_type;
	using splitter_type = hyper_plane_type;
	using splitter_seq_type = hyper_plane_seq_type;
	// using AugType = std::optional<bool>;
	using split_rule_type = SplitRule;

	// NOTE: for the CPAM
	using curve_code_type = typename Point::at_type::curve_code_type;
	using id_type = typename Point::at_type::id_type;

	using cpam_key_type = typename Point::at_type; // morton_id, id
	using cpam_val_type = coords_type;
	// using cpam_aug_type = std::pair<box_type, size_t>;
	// PERF: aug value can remove the size, because the tree has this info
	using cpam_aug_type = box_type;

	struct cpam_entry {
		using key_t = cpam_key_type;
		using val_t = cpam_val_type;
		using aug_t = cpam_aug_type;
		using entry_t = Point;

		using filling_curve_t = SplitRule;
		using key_entry_pointer =
			std::pair<typename Point::at_type, entry_t *>;
		// using entry_t_ref_v = Point*;
		using entry_t_ref_wrapper_v = std::reference_wrapper<Point>;

		static inline key_t const &get_key(entry_t const &e)
		{
			return e.get_aug();
		}
		static inline val_t const &get_val(entry_t const &e)
		{
			return e.get_coords();
		}
		static inline void set_val(entry_t &e, val_t const &v)
		{
			e.get_coords() = v;
		}
		static inline entry_t to_entry(key_t const &k, val_t const &v)
		{
			return entry_t(v, k);
			// return std::make_tuple(k, v);
		};
		static inline aug_t from_entry(entry_t const &e)
		{
			return from_entry(get_key(e), get_val(e));
		}

		// how to compare key
		static inline bool comp(key_t const &a, key_t const &b)
		{
			return a < b;
		}

		// get an empty aug
		static aug_t get_empty()
		{
			return base_type::get_empty_box();
		}

		// WARN: this invoke implicity conversion from coords_type to
		// basic_point
		static aug_t from_entry(key_t const &, val_t const &v)
		{
			return box_type(v, v);
		}

		static aug_t from_entry_array(entry_t *et, size_t sz)
		{
			return base_type::get_box(parlay::slice(et, et + sz));
		}

		// combine two aug val to get a new aug
		static aug_t combine(aug_t const &a, aug_t const &b)
		{
			return base_type::get_box(a, b);
		}
	};

	// using cpam_aug_map_type = cpam::aug_map<cpam_entry,
	// base_type::leaf_capacity>;
	using cpam_aug_map_type = cpam::aug_map<cpam_entry, 40>;
	using cpam_map_type = cpam_aug_map_type::map_type;

	using cpam_inner_entry_type =
		std::tuple<typename cpam_entry::key_t,
			   std::reference_wrapper<typename cpam_entry::val_t>>;

	using leaf_type = cpam_aug_map_type::Tree::node;
	using interior_type = cpam_aug_map_type::Tree::node;

	// NOTE: general tree structure
	/* delete_tree_wrapper is idempotent, so an explicit delete_tree()
	 * before this stays correct. */
	~p_tree() override
	{
		delete_tree();
	}

	p_tree() = default;
	/*
	 * Move-only: root_ is an owning raw pointer. Assignment ends up
	 * deleted because SplitRule holds const members; construction is
	 * what containers and factory returns need.
	 */
	p_tree(p_tree &&) = default;
	p_tree &operator=(p_tree &&) = default;

	void p_tree_tag();

	template <typename Range>
	void build(Range &&in);

	template <typename Range>
	void batch_insert(Range &&in);

	// NOTE: every point is assumed to be in the tree; if that may not
	// hold, use batch_diff
	template <typename Range>
	void batch_delete(Range &&in);

	// NOTE: tolerates points that are not in the tree
	template <typename Range>
	void batch_diff(Range &&in);

	/*
	 * Copies every point into out, which must be exactly get_size()
	 * long: the recursion splits out by subtree size, so any other
	 * length would write outside it. Returns the number written, or 0
	 * if out was the wrong size and nothing was written.
	 */
	template <typename Range>
	size_t flatten(Range &&out) const;

	template <typename Range>
	/* Not const: cpam takes the map by non-const reference. */
	auto knn(Point const &q, bounded_queue<Point, Range> &bq);

	auto range_count(box_type const &query_box);

	template <typename Range>
	auto range_query(box_type const &query_box, Range &&out);

	constexpr void delete_tree() override;

	size_t get_size() const
	{
		return cpam_aug_map_.size();
	}

	constexpr static char const *get_tree_name()
	{
		return "PTree";
	}

	constexpr static char const *check_has_box()
	{
		return "HasBox";
	}

private:
	void build_(slice_type in);
	void batch_insert_(slice_type in);
	void batch_delete_(slice_type in);
	void batch_diff_(slice_type in);

	SplitRule space_filling_curve_;
	cpam_aug_map_type cpam_aug_map_;
};

} // namespace psi

#include "p_tree_impl/p_batch_delete.hpp"
#include "p_tree_impl/p_batch_diff.hpp"
#include "p_tree_impl/p_batch_insert.hpp"
#include "p_tree_impl/p_build_tree.hpp"
#include "p_tree_impl/p_inter_node.hpp"
#include "p_tree_impl/p_override.hpp"

#endif // PSI_P_TREE_H
