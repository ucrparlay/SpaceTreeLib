#ifndef PSI_ORTH_TREE_H
#define PSI_ORTH_TREE_H

#include <functional>
#include <utility>

#include "base_tree.h"
#include "dependence/augmentation.h"

namespace psi
{

template <typename Point, typename SplitRule,
	  typename LeafAugType = box_leaf_aug<base_tree<Point>>,
	  typename InteriorAugType = box_interior_aug<base_tree<Point>>,
	  uint_fast8_t md = 2, uint_fast8_t SkHeight = 6,
	  uint_fast8_t ImbaRatio = 30>
class orth_tree
    : public base_tree<Point,
		       orth_tree<Point, SplitRule, LeafAugType, InteriorAugType,
				 md, SkHeight, ImbaRatio>,
		       SkHeight, ImbaRatio>
{
public:
	static constexpr size_t splitter_num = md;
	static constexpr size_t node_regions = 1 << md;

	using base_type =
		base_tree<Point,
			  orth_tree<Point, SplitRule, LeafAugType,
				    InteriorAugType, md, SkHeight, ImbaRatio>,
			  SkHeight, ImbaRatio>;

	static_assert(
		leaf_augmentation<LeafAugType, typename base_type::slice_type>,
		"LeafAugType needs A(), A(slice_type), update_aug(slice_type), "
		"reset()");
	static_assert(interior_augmentation<InteriorAugType>,
		      "InteriorAugType needs set_parallel_flag(bool), "
		      "reset_parallel_flag(), get_parallel_flag_ini_status(), "
		      "force_parallel(size_t)");
	/*
	 * divide_rotate hands split_sample a null slice, which only
	 * spatial_median ignores. object_median dereferences it and the build
	 * segfaults.
	 */
	static_assert(is_spatial_median_split<
			      typename SplitRule::partition_rule_type>,
		      "orth_tree needs a spatial_median partition rule");

	using bucket_type = base_type::bucket_type;
	using balls_type = base_type::balls_type;
	using bucket_seq_type = base_type::bucket_seq_type;
	using dims_type = base_type::dims_type;
	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	using num_type = num_comparator<coord_type>;
	using slice_type = base_type::slice_type;
	using points_type = base_type::points_type;
	using points_iter_type = base_type::points_iter_type;
	using box_type = base_type::box_type;
	using box_seq_type = base_type::box_seq_type;
	// using circle_type = base_type::circle_type;

	using hyper_plane_type = base_type::hyper_plane_type;
	using hyper_plane_seq_type = base_type::hyper_plane_seq_type;
	using splitter_type = std::array<hyper_plane_type, splitter_num>;
	using splitter_seq_type = parlay::sequence<splitter_type>;
	using node_box_seq_type = base_type::node_box_seq_type;
	using node_box_type = base_type::node_box_type;
	// using AugType = std::optional<bool>;

	// struct kd_interior_node;
	struct orth_interior_node;

	using split_rule_type = SplitRule;
	using leaf_type = leaf_node<Point, slice_type, base_type::leaf_capacity,
				    LeafAugType, parlay::move_assign_tag>;
	using interior_type = orth_interior_node;
	using orth_node_arr_type = interior_type::orth_node_arr_type;
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
	friend SplitRule;

	/* delete_tree_wrapper is idempotent, so an explicit delete_tree()
	 * before this stays correct. */
	~orth_tree() override
	{
		delete_tree();
	}

	void orth_tree_tag();

	template <typename Range, typename... Args>
	void build(Range &&in, Args &&...args);

	template <typename Range>
	void batch_insert(Range &&in);

	// NOTE: every point is assumed to be in the tree; if that may not
	// hold, use batch_diff
	template <typename Range>
	void batch_delete(Range &&in);

	// NOTE: tolerates points that are not in the tree
	template <typename Range>
	void batch_diff(Range &&in);

	template <typename Range>
	void flatten(Range &&out) const;

	template <typename Range>
	auto knn(Point const &q, bounded_queue<Point, Range> &bq) const;

	auto range_count(box_type const &query_box) const;

	template <typename Range>
	auto range_query(box_type const &query_box, Range &&out) const;

	constexpr void delete_tree() override;

	void set_bounding_box(box_type const &box)
	{
		this->tree_box_ = box;
		fixed_box = true;
	}

	constexpr static char const *get_tree_name()
	{
		return "OrthTree";
	}

	constexpr static char const *check_has_box()
	{
		if constexpr (has_box<InteriorAugType>)
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

	SplitRule split_rule_;
	bool fixed_box = false;
};

} // namespace psi

#include "orth_tree_impl/orth_batch_delete.hpp"
#include "orth_tree_impl/orth_batch_diff.hpp"
#include "orth_tree_impl/orth_batch_insert.hpp"
#include "orth_tree_impl/orth_build_tree.hpp"
#include "orth_tree_impl/orth_inter_node.hpp"
#include "orth_tree_impl/orth_override.hpp"

#endif // PSI_ORTH_TREE_H
