#ifndef PSI_KD_TREE_H
#define PSI_KD_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "psi/base_tree.h"
#include "psi/dependence/augmentation.h"
#include "psi/dependence/concepts.h"

namespace psi
{

template <typename Point, typename SplitRule,
	  typename LeafAugType = box_leaf_aug<base_tree<Point>>,
	  typename InteriorAugType = box_interior_aug<base_tree<Point>>,
	  uint_fast8_t SkHeight = 6, uint_fast8_t ImbaRatio = 30>
class kd_tree : public base_tree<Point,
				 kd_tree<Point, SplitRule, LeafAugType,
					 InteriorAugType, SkHeight, ImbaRatio>,
				 SkHeight, ImbaRatio>
{
public:
	using base_type =
		base_tree<Point,
			  kd_tree<Point, SplitRule, LeafAugType,
				  InteriorAugType, SkHeight, ImbaRatio>,
			  SkHeight, ImbaRatio>;

	static_assert(
		leaf_augmentation<LeafAugType, typename base_type::slice_type>,
		"LeafAugType needs A(), A(slice_type), update_aug(slice_type), "
		"reset()");
	static_assert(interior_augmentation<InteriorAugType>,
		      "InteriorAugType needs set_parallel_flag(bool), "
		      "reset_parallel_flag(), get_parallel_flag_ini_status(), "
		      "force_parallel(size_t)");

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
	using split_rule_type = SplitRule;

	static constexpr dims_type const md = 2;

	using leaf_type = leaf_node<Point, slice_type, base_type::leaf_capacity,
				    LeafAugType, parlay::move_assign_tag>;
	struct kd_interior_node;
	using interior_type = kd_interior_node;

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
	~kd_tree() override
	{
		delete_tree();
	}

	kd_tree() = default;
	/*
	 * Move-only: root_ is an owning raw pointer. Assignment ends up
	 * deleted because SplitRule holds const members; construction is
	 * what containers and factory returns need.
	 */
	kd_tree(kd_tree &&) = default;
	kd_tree &operator=(kd_tree &&) = default;

	void kd_tree_tag();

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
	auto knn(Point const &q, bounded_queue<Point, Range> &bq) const;

	auto range_count(box_type const &query_box) const;

	template <typename Range>
	auto range_query(box_type const &query_box, Range &&out) const;

	constexpr void delete_tree() override;

	constexpr static char const *get_tree_name()
	{
		return "KdTree";
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
	node *build_recursive(slice_type in, slice_type out, dims_type dim,
			      box_type const &bx);
	node *serial_build_recursive(slice_type in, slice_type out,
				     dims_type dim, box_type const &bx);
	void divide_rotate(slice_type in, splitter_seq_type &pivots,
			   dims_type dim, bucket_type idx,
			   box_seq_type &box_seq, box_type const &bx);
	void pick_pivots(slice_type in, size_t const &n,
			 splitter_seq_type &pivots, dims_type const dim,
			 box_seq_type &box_seq, box_type const &bx);

	void batch_insert_(slice_type in);
	node *batch_insert_recursive(node *T, slice_type in, slice_type out,
				     dims_type d);

	void batch_delete_(slice_type in);
	node_box_type batch_delete_recursive(node *T, box_type const &bx,
					     slice_type in, slice_type out,
					     dims_type d, bool has_tomb);

	void batch_diff_(slice_type in);
	node_box_type batch_diff_recursive(node *T, box_type const &bx,
					   slice_type in, slice_type out,
					   dims_type d);

	SplitRule split_rule_;
};

} // namespace psi

#include "psi/kd_tree_impl/kd_batch_delete.hpp"
#include "psi/kd_tree_impl/kd_batch_diff.hpp"
#include "psi/kd_tree_impl/kd_batch_insert.hpp"
#include "psi/kd_tree_impl/kd_build_tree.hpp"
#include "psi/kd_tree_impl/kd_inter_node.hpp"
#include "psi/kd_tree_impl/kd_override.hpp"

#endif // PSI_KD_TREE_H
