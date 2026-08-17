#ifndef PSI_P_TREE_H
#define PSI_P_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "psi/base_tree.h"
#include "psi/tree_traits.h"

/*
 * cpam is third party. Silence its noise, but not -Wreturn-local-addr or
 * -Wuninitialized: both currently fire in augmented_node.h and both are real.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "psi/dependence/cpam/cpam.h"
#pragma GCC diagnostic pop

namespace psi
{

template <typename Traits>
class p_tree : public base_tree<Traits, p_tree<Traits>>
{
public:
	using base_type = base_tree<Traits, p_tree<Traits>>;

	using point_type = typename base_type::point_type;
	using bucket_type = typename base_type::bucket_type;
	using balls_type = typename base_type::balls_type;
	using dims_type = typename base_type::dims_type;
	using bucket_seq_type = typename base_type::bucket_seq_type;
	using coord_type = typename base_type::coord_type;
	using coords_type = typename base_type::coords_type;
	using num_type = typename base_type::num_type;
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
	using split_rule_type = typename base_type::split_rule_type;

	/* for the CPAM */
	using curve_code_type = typename point_type::at_type::curve_code_type;
	using id_type = typename point_type::at_type::id_type;

	using cpam_key_type = typename point_type::at_type; /* morton_id, id */
	using cpam_val_type = coords_type;

	using cpam_aug_type = box_type;

	struct cpam_entry {
		using key_t = cpam_key_type;
		using val_t = cpam_val_type;
		using aug_t = cpam_aug_type;
		using entry_t = point_type;

		using filling_curve_t = split_rule_type;
		using key_entry_pointer =
			std::pair<typename point_type::at_type, entry_t *>;
		using entry_t_ref_wrapper_v =
			std::reference_wrapper<point_type>;

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
		};
		static inline aug_t from_entry(entry_t const &e)
		{
			return from_entry(get_key(e), get_val(e));
		}

		/* how to compare key */
		static inline bool comp(key_t const &a, key_t const &b)
		{
			return a < b;
		}

		/* get an empty aug */
		static aug_t get_empty()
		{
			return base_type::get_empty_box();
		}

		/*
		 * this invoke implicity conversion from coords_type to
		 * basic_point
		 */
		static aug_t from_entry(key_t const &, val_t const &v)
		{
			return box_type(v, v);
		}

		static aug_t from_entry_array(entry_t *et, size_t sz)
		{
			return base_type::get_box(parlay::slice(et, et + sz));
		}

		/* combine two aug val to get a new aug */
		static aug_t combine(aug_t const &a, aug_t const &b)
		{
			return base_type::get_box(a, b);
		}
	};

	using cpam_aug_map_type = cpam::aug_map<cpam_entry, 40>;
	using cpam_map_type = typename cpam_aug_map_type::map_type;

	using cpam_inner_entry_type =
		std::tuple<typename cpam_entry::key_t,
			   std::reference_wrapper<typename cpam_entry::val_t>>;

	using leaf_type = typename cpam_aug_map_type::Tree::node;
	using interior_type = typename cpam_aug_map_type::Tree::node;

	/* general tree structure */
	/* delete_tree_wrapper is idempotent, so an explicit delete_tree()
	 * before this stays correct. */
	~p_tree() override
	{
		delete_tree();
	}

	p_tree() = default;
	/*
	 * Move-only: root_ is an owning raw pointer. Assignment ends up
	 * deleted because the split rule holds const members; construction is
	 * what containers and factory returns need.
	 */
	p_tree(p_tree &&) = default;
	p_tree &operator=(p_tree &&) = default;

	void p_tree_tag();

	template <typename Range>
	void build(Range &&in);

	template <typename Range>
	void batch_insert(Range &&in);

	/*
	 * every point is assumed to be in the tree; if that may not
	 * hold, use batch_diff
	 */
	template <typename Range>
	void batch_delete(Range &&in);

	/* tolerates points that are not in the tree */
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
	auto knn(point_type const &q,
		 bounded_queue<point_type, Range> &bq) const;

	/* The buffer-taking forms below would otherwise hide the
	 * allocating ones the base provides. */
	using base_type::empty;
	using base_type::knn;
	using base_type::range_query;

	auto range_count(box_type const &query_box) const;

	template <typename Range>
	auto range_query(box_type const &query_box, Range &&out) const;

	constexpr void delete_tree() override;

	size_t get_size() const
	{
		return cpam_aug_map_.size();
	}

	/*
	 * Not the base's: p_tree never fills tree_box_. Its bounding box is
	 * the cpam augmentation at the root, which the map keeps current
	 * through every update, so unlike the base's this one is always tight.
	 */
	box_type bounds() const
	{
		return cpam_aug_map_.aug_val();
	}

	/*
	 * Parsed by script_ae/merge_*.py to label the figures; changing it
	 * silently relabels published results.
	 */
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

	/*
	 * No split-rule member: unlike kd_tree and orth_tree, p_tree never
	 * calls its rule. The rule reaches the tree as cpam_entry's
	 * filling_curve_t, and the code is applied by cpam_sample_sort during
	 * the build. Adding a curve means writing one with encode(), not
	 * touching anything here.
	 *
	 * The queries do not change the tree, but cpam takes its map by
	 * non-const reference.
	 */
	mutable cpam_aug_map_type cpam_aug_map_;
};

} /* namespace psi */

#include "psi/p_tree_impl/p_batch_delete.hpp"
#include "psi/p_tree_impl/p_batch_diff.hpp"
#include "psi/p_tree_impl/p_batch_insert.hpp"
#include "psi/p_tree_impl/p_build_tree.hpp"
#include "psi/p_tree_impl/p_inter_node.hpp"
#include "psi/p_tree_impl/p_override.hpp"

#endif /* PSI_P_TREE_H */
