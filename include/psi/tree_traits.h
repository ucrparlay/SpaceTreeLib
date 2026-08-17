#ifndef PSI_TREE_TRAITS_H_
#define PSI_TREE_TRAITS_H_

#include <cstdint>
#include <tuple>
#include <type_traits>
#include <utility>

#include "psi/dependence/augmentation.h"
#include "psi/dependence/comparator.h"
#include "psi/dependence/concepts.h"
#include "psi/dependence/tree_node.h"

namespace psi
{

/*
 * The vocabulary that follows from the point type alone: every alias a tree
 * needs, the constants that tune its shape, and the box arithmetic all three
 * trees share.
 *
 * It is a separate type from tree_traits because the default augmentations
 * need the box operations, and they are named in tree_traits' own parameter
 * list -- a template cannot refer to itself there.
 */
template <typename Point>
struct point_traits {
	using point_type = Point;
	using basic_point = typename Point::bp_type;

	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	using dis_type = typename Point::dis_type;
	using num_type = num_comparator<coord_type>;

	using dims_type = uint_fast8_t;
	using balls_type = uint_fast32_t;
	using depth_type = int;
	using ball_seq_type = parlay::sequence<balls_type>;

	using slice_type = parlay::slice<Point *, Point *>;
	using points_type = parlay::sequence<Point>;
	using points_iter_type = typename parlay::sequence<Point>::iterator;

	using hyper_plane_type = std::pair<coord_type, dims_type>;
	using hyper_plane_seq_type = parlay::sequence<hyper_plane_type>;
	using box_type = std::pair<basic_point, basic_point>;
	using box_seq_type = parlay::sequence<box_type>;

	/* What the convenience knn hands back: points by value, so the result
	 * outlives any update to the tree. */
	using knn_result_type = std::pair<Point, dis_type>;
	using knn_result_seq_type = parlay::sequence<knn_result_type>;

	using node_box_type = std::pair<node *, box_type>;
	using node_box_seq_type = parlay::sequence<node_box_type>;
	using node_tag_type = std::pair<node *, uint_fast8_t>;
	using node_tag_seq_type = parlay::sequence<node_tag_type>;

	static constexpr dims_type const num_dims =
		std::tuple_size_v<coords_type>;

	/* tree structure */
	static constexpr uint_fast8_t const leaf_capacity = 32;
	static constexpr uint_fast8_t const sparse_leaf_threshold = 24;
	static constexpr uint_fast16_t const serial_build_cutoff = 1 << 10;

	/* block param in partition */
	static constexpr uint_fast8_t const log2_base = 10;
	static constexpr uint_fast16_t const block_size = 1 << log2_base;

	/* box_type operations */
	static inline coord_type get_box_mid(dims_type const d,
					     box_type const &bx);
	static inline bool legal_box(box_type const &bx);
	/*
	 * True when the box covers nothing, i.e. some dimension has first >
	 * second. Intersecting two boxes produces one, so it is a legal query
	 * whose answer is zero -- the public queries return early rather than
	 * asserting their way down the traversal.
	 */
	static inline bool box_is_empty(box_type const &bx);
	static inline bool within_box(box_type const &a, box_type const &b);
	static inline bool same_box(box_type const &a, box_type const &b);
	static inline bool within_box(Point const &p, box_type const &bx);
	static inline bool box_intersect_box(box_type const &a,
					     box_type const &b);
	static inline bool is_box_line_in_dimension(box_type const &box,
						    dims_type d);
	static inline bool vertical_line_split_box(coord_type const &l,
						   box_type const &box,
						   dims_type d);
	static inline bool vertical_line_on_box_left_edge(coord_type const &l,
							  box_type const &box,
							  dims_type d);
	static inline bool vertical_line_on_box_right_edge(coord_type const &l,
							   box_type const &box,
							   dims_type d);
	static inline bool vertical_line_on_box_edge(coord_type const &l,
						     box_type const &box,
						     dims_type d);
	static inline bool vertical_line_intersect_box(coord_type const &l,
						       box_type const &box,
						       dims_type d);
	static inline bool
	vertical_line_intersect_box_exclude(coord_type const &l,
					    box_type const &box, dims_type d);

	static inline box_type get_empty_box();
	static inline Point get_box_center(box_type const &box);
	static box_type get_box(box_type const &x, box_type const &y);
	static box_type get_box(slice_type V);
	static box_type get_box(box_seq_type const &box_seq);
	template <typename leaf_type, typename interior_type>
	static box_type get_box(node *T);
};

/*
 * Everything a tree is configured with, in one type. Name it once and hand
 * it to any of the trees:
 *
 *     using traits = psi::tree_traits<point, split_rule>;
 *     psi::kd_tree<traits> kd;
 *     psi::orth_tree<psi::orth_tree_traits<point, split_rule>> orth;
 */
template <typename Point, typename SplitRule,
	  typename LeafAug = box_leaf_aug<point_traits<Point>>,
	  typename InteriorAug = box_interior_aug<point_traits<Point>>,
	  uint_fast8_t SkHeight = 6, uint_fast8_t ImbaRatio = 30>
struct tree_traits : point_traits<Point> {
	using split_rule_type = SplitRule;
	using leaf_aug_type = LeafAug;
	using interior_aug_type = InteriorAug;

	static_assert(
		leaf_augmentation<LeafAug,
				  typename point_traits<Point>::slice_type>,
		"LeafAug needs A(), A(slice_type), update_aug(slice_type), "
		"reset()");
	static_assert(interior_augmentation<InteriorAug>,
		      "InteriorAug needs set_parallel_flag(bool), "
		      "reset_parallel_flag(), get_parallel_flag_ini_status(), "
		      "force_parallel(size_t)");

	/*
	 * when SkHeight >= 8, the # bucket is 255, total skeleton nodes
	 * >= 255*2
	 */
	using bucket_type =
		std::conditional_t<(SkHeight > 7), uint_fast16_t, uint_fast8_t>;
	using bucket_seq_type = parlay::sequence<bucket_type>;

	/*
	 * uint32t handle up to 4e9 at least
	 * bucket num should smaller than 1<<8 to handle type overflow
	 */
	static constexpr bucket_type const build_depth_once = SkHeight;
	static constexpr bucket_type const pivot_num =
		(1 << build_depth_once) - 1;
	static constexpr bucket_type const bucket_num = 1 << build_depth_once;

	/* reconstruct weight threshold */
	static constexpr uint_fast8_t const imbalance_ratio = ImbaRatio;
};

/*
 * Skeleton height for a given dimension. A skeleton level splits every
 * dimension once, so the height has to be a multiple of the dimension.
 * Throwing is how a consteval function reports an unsupported dimension.
 */
consteval uint_fast8_t orth_build_depth_once(uint_fast8_t const dim)
{
	if (dim == 2 || dim == 3) {
		return 6;
	}
	if (dim == 4) {
		return 8;
	}
	if (dim >= 5 && dim <= 8) {
		return dim;
	}
	throw "orth_tree has no skeleton height for this dimension";
}

/* tree_traits with the skeleton height orth_tree needs for the dimension. */
template <typename Point, typename SplitRule,
	  typename LeafAug = box_leaf_aug<point_traits<Point>>,
	  typename InteriorAug = box_interior_aug<point_traits<Point>>,
	  uint_fast8_t SkHeight = orth_build_depth_once(Point::get_dim()),
	  uint_fast8_t ImbaRatio = 30>
using orth_tree_traits = tree_traits<Point, SplitRule, LeafAug, InteriorAug,
				     SkHeight, ImbaRatio>;

} /* namespace psi */

#include "psi/tree_traits_impl/box_op.hpp"

#endif /* PSI_TREE_TRAITS_H_ */
