#ifndef PSI_ORTH_TREE_H
#define PSI_ORTH_TREE_H

#include <functional>
#include <utility>

#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"
#include "psi/tree_traits.h"

namespace psi
{

template <typename Traits>
class orth_tree : public base_tree<Traits, orth_tree<Traits>>
{
public:
	using base_type = base_tree<Traits, orth_tree<Traits>>;

	using point_type = typename base_type::point_type;
	using bucket_type = typename base_type::bucket_type;
	using balls_type = typename base_type::balls_type;
	using bucket_seq_type = typename base_type::bucket_seq_type;
	using dims_type = typename base_type::dims_type;
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
	using node_box_seq_type = typename base_type::node_box_seq_type;
	using node_box_type = typename base_type::node_box_type;
	using split_rule_type = typename base_type::split_rule_type;

	/* One split per dimension, so a node has 2^dim children. */
	static constexpr dims_type const md = base_type::num_dims;
	static constexpr size_t splitter_num = md;
	static constexpr size_t node_regions = 1 << md;

	using splitter_type = std::array<hyper_plane_type, splitter_num>;
	using splitter_seq_type = parlay::sequence<splitter_type>;

	/*
	 * divide_rotate hands split_sample a null slice, which only
	 * spatial_median ignores. object_median dereferences it and the build
	 * segfaults.
	 */
	static_assert(is_spatial_median_split<
			      typename split_rule_type::partition_rule_type>,
		      "orth_tree needs a spatial_median partition rule");
	/*
	 * A skeleton level splits every dimension once, so a height that is
	 * not a multiple of the dimension leaves the skeleton lopsided.
	 * orth_tree_traits picks a height that fits.
	 */
	static_assert(base_type::build_depth_once % md == 0,
		      "orth_tree's skeleton height must be a multiple of the "
		      "point's dimension -- use psi::orth_tree_traits");

	struct orth_interior_node;

	using leaf_type =
		leaf_node<point_type, slice_type, base_type::leaf_capacity,
			  typename Traits::leaf_aug_type,
			  parlay::move_assign_tag>;
	using interior_type = orth_interior_node;
	using node_arr_type = typename interior_type::node_arr_type;
	using inner_tree =
		typename base_type::template inner_tree<leaf_type,
							interior_type>;
	using box_cut_type = typename base_type::box_cut_type;

	template <typename leaf_type, typename interior_type, bool granularity,
		  typename... Args>
	friend node *base_type::rebuild_single_tree(node *T, Args &&...args);

	template <typename leaf_type, typename interior_type,
		  typename PrepareFunc, typename... Args>
	friend node *
	base_type::rebuild_with_insert(node *T, PrepareFunc prepare_func,
				       slice_type in, Args &&...args);

	/* The split rule calls build_recursive back; see divide_space. */
	friend split_rule_type;

	/* delete_tree_wrapper is idempotent, so an explicit delete_tree()
	 * before this stays correct. */
	~orth_tree() override
	{
		delete_tree();
	}

	orth_tree() = default;
	/*
	 * Move-only: root_ is an owning raw pointer. Assignment ends up
	 * deleted because the split rule holds const members; construction is
	 * what containers and factory returns need.
	 */
	orth_tree(orth_tree &&) = default;
	orth_tree &operator=(orth_tree &&) = default;

	void orth_tree_tag();

	template <typename Range, typename... Args>
	void build(Range &&in, Args &&...args);

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

	void set_bounding_box(box_type const &box)
	{
		this->tree_box_ = box;
		fixed_box = true;
	}

	/*
	 * Parsed by script_ae/merge_*.py to label the figures; changing it
	 * silently relabels published results.
	 */
	constexpr static char const *get_tree_name()
	{
		return "OrthTree";
	}

	constexpr static char const *check_has_box()
	{
		if constexpr (has_box<typename Traits::interior_aug_type>)
			return "HasBox";
		else
			return "NoBox";
	}

private:
	void build_(slice_type in);
	void build_(slice_type in, box_type const &box);
	node *build_recursive(slice_type in, slice_type out,
			      box_type const &bx);
	node *serial_build_recursive(slice_type in, slice_type out,
				     box_type const &bx,
				     bool checked_duplicate);
	void serial_split(slice_type in, dims_type dim, dims_type idx,
			  box_type const &box,
			  parlay::sequence<balls_type> &sums);
	void serial_split_skeleton(node *T, slice_type in, dims_type dim,
				   dims_type idx,
				   parlay::sequence<balls_type> &sums);
	void divide_rotate(hyper_plane_seq_type &pivots, dims_type dim,
			   bucket_type idx, box_seq_type &box_seq,
			   box_type const &box);
	void pick_pivots(slice_type in, size_t const &n,
			 hyper_plane_seq_type &pivots, dims_type const dim,
			 box_seq_type &box_seq, box_type const &bx);

	void batch_insert_(slice_type in);
	node *batch_insert_recursive(node *T, slice_type in, slice_type out,
				     box_type const &bx);

	void batch_delete_(slice_type in);
	node *batch_delete_recursive(node *T, slice_type in, slice_type out,
				     box_type const &box, bool has_tomb);

	void batch_diff_(slice_type in);
	node *batch_diff_recursive(node *T, slice_type in, slice_type out);

	split_rule_type split_rule_;
	bool fixed_box = false;
};

} /* namespace psi */

#include "psi/orth_tree_impl/orth_batch_delete.hpp"
#include "psi/orth_tree_impl/orth_batch_diff.hpp"
#include "psi/orth_tree_impl/orth_batch_insert.hpp"
#include "psi/orth_tree_impl/orth_build_tree.hpp"
#include "psi/orth_tree_impl/orth_inter_node.hpp"
#include "psi/orth_tree_impl/orth_override.hpp"

#endif /* PSI_ORTH_TREE_H */
