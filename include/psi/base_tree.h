#ifndef PSI_BASE_TREE_H_
#define PSI_BASE_TREE_H_

#include <cstdint>
#include <cstdio>
#include <type_traits>
#include <utility>

#include "psi/dependence/comparator.h"
#include "psi/dependence/loggers.h"
#include "psi/dependence/search_container.h"
#include "psi/dependence/tree_node.h"

namespace psi
{

// NOTE: Basetree
template <typename Point, typename DerivedTree = void,
	  uint_fast8_t SkHeight = 6, uint_fast8_t ImbaRatio = 30>
class base_tree
{
public:
	// NOTE: when SkHeight >= 8, the # bucket is 255, total skeleton nodes
	// >= 255*2

	using basic_point = Point::bp_type;
	using template_point_type = Point;

	using bucket_type =
		std::conditional_t<(SkHeight > 7), uint_fast16_t, uint_fast8_t>;
	using balls_type = uint_fast32_t;
	using dims_type = uint_fast8_t;
	// using depth_type = uint_fast8_t;
	// using depth_type = int_fast8_t;
	using depth_type = int;
	using bucket_seq_type = parlay::sequence<bucket_type>;
	using ball_seq_type = parlay::sequence<balls_type>;

	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	using dis_type = typename Point::dis_type;
	using num_type = num_comparator<coord_type>;
	using slice_type = parlay::slice<Point *, Point *>;
	using points_type = parlay::sequence<Point>;
	using points_iter_type = typename parlay::sequence<Point>::iterator;
	using hyper_plane_type = std::pair<coord_type, dims_type>;
	using hyper_plane_seq_type = parlay::sequence<hyper_plane_type>;
	using box_type = std::pair<basic_point, basic_point>;
	using box_seq_type = parlay::sequence<box_type>;

	using node_box_type = std::pair<node *, box_type>;
	using node_box_seq_type = parlay::sequence<node_box_type>;
	using node_tag_type = std::pair<node *, uint_fast8_t>;
	using node_tag_seq_type = parlay::sequence<node_tag_type>;

	// TODO: use the one provided in aug

	// NOTE: Const variables
	// NOTE: uint32t handle up to 4e9 at least
	// WARN: bucket num should smaller than 1<<8 to handle type overflow
	static constexpr dims_type const num_dims =
		std::tuple_size_v<coords_type>;
	static constexpr bucket_type const build_depth_once = SkHeight;
	static constexpr bucket_type const pivot_num =
		(1 << build_depth_once) - 1;
	static constexpr bucket_type const bucket_num = 1 << build_depth_once;

	// NOTE: tree structure
	static constexpr uint_fast8_t const leaf_capacity = 32;
	static constexpr uint_fast8_t const sparse_leaf_threshold = 24;
	// static constexpr uint_fast8_t const leaf_capacity = 8;
	// static constexpr uint_fast8_t const sparse_leaf_threshold = 4;
	static constexpr uint_fast16_t const serial_build_cutoff = 1 << 10;

	// NOTE: block param in partition
	static constexpr uint_fast8_t const log2_base = 10;
	static constexpr uint_fast16_t const block_size = 1 << log2_base;

	// NOTE: reconstruct weight threshold
	static constexpr uint_fast8_t const imbalance_ratio = ImbaRatio;

	// NOTE: array based inner tree for batch insertion and deletion
	template <typename leaf_type, typename interior_type>
	struct inner_tree;

	// NOTE: compute the bounding box on the fly
	struct box_cut_type;

	// NOTE: get the imbalance ratio
	static inline size_t get_imbalance_ratio();
	static inline bool imbalance_node(size_t const l, size_t const n);
	static inline bool sparse_node(size_t const l, size_t const n);

	// NOTE: box_type operations
	static inline coord_type get_box_mid(dims_type const d,
					     box_type const &bx);
	static inline bool legal_box(box_type const &bx);
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
	template <typename leaf_type, typename interior_type>
	static inline decltype(auto) retrieve_box(node const *T)
		requires(has_box<typename leaf_type::at_type> &&
			 has_box<typename interior_type::at_type>);

	static inline box_type get_empty_box();
	static inline Point get_box_center(box_type const &box);
	static box_type get_box(box_type const &x, box_type const &y);
	static box_type get_box(slice_type V);
	static box_type get_box(box_seq_type const &box_seq);
	template <typename leaf_type, typename interior_type>
	static box_type get_box(node *T);

	// NOTE: circle_type operations
	struct normal_circle {
		Point const &get_center() const
		{
			return center;
		}

		coord_type get_radius() const
		{
			return radius;
		}

		coord_type get_radius_square() const
		{
			return radius * radius;
		}

		static coord_type compute_radius(double const &r)
		{
			return static_cast<coord_type>(r);
		}

		friend std::ostream &operator<<(std::ostream &o,
						normal_circle const &cl)
		{
			o << "{ " << cl.center << ", " << cl.radius << "}";
			return o;
		}

		Point center;
		coord_type radius;
	};

	struct cover_circle {
		Point const &get_center() const
		{
			return center;
		}

		coord_type get_radius() const
		{
			return static_cast<coord_type>(1 << level);
		}

		coord_type get_radius_square() const
		{
			return get_radius() * get_radius();
		}

		static depth_type compute_radius(auto const &r)
		{
			// TODO: hard to represent the radius with length 0
			return num_comparator<decltype(r)>::is_zero(r)
				       ? -1
				       : static_cast<depth_type>(
						 std::ceil(std::log2(r)));
		}

		bool operator==(cover_circle const &cl) const
		{
			return center == cl.center && level == cl.level;
		}

		friend std::ostream &operator<<(std::ostream &o,
						cover_circle const &cl)
		{
			o << "{ " << cl.center << ", " << cl.level << "}";
			return o;
		}

		Point center;
		depth_type level;
	};

	template <typename CircleType>
	static bool legal_circle(CircleType const &cl);

	template <typename CircleType>
	static inline bool within_circle(Point const &p, CircleType const &cl);

	template <typename CircleType>
	static inline bool within_circle(box_type const &box,
					 CircleType const &cl);

	template <typename CircleType1, typename CircleType2>
	static inline bool circle_within_circle(CircleType1 const &a,
						CircleType2 const &b);

	template <typename CircleType>
	static inline bool circle_intersect_box(CircleType const &cl,
						box_type const &box);

	template <typename CircleType>
	static inline bool circle_intersect_circle(CircleType const &a,
						   CircleType const &b);

	template <typename CircleType>
	static inline CircleType get_circle(box_type const &box);

	template <typename CircleType>
	static inline CircleType get_circle(slice_type V);

	template <typename CircleType>
	static inline CircleType get_circle(CircleType const &a,
					    CircleType const &b);

	template <typename CircleType>
	static inline CircleType get_circle(Point const &p,
					    CircleType const &cl);

	template <typename CircleType>
	static inline CircleType
	get_circle(parlay::sequence<CircleType> const &circle_seq);

	/*
	 * Copy a caller range into tree owned storage and hand the slice to op.
	 * Every public build, batch_insert, batch_delete and batch_diff funnels
	 * through here. The copy is unconditional: the trees consume the buffer
	 * in place, and a move-when-rvalue path would change how much memory
	 * traffic every measured operation does.
	 *
	 * Preconditions that read the caller's own range must be checked by the
	 * caller, before this runs.
	 */
	template <typename Range, typename Fn>
	static void ingest_range(Range &&in, Fn &&op)
	{
		static_assert(parlay::is_random_access_range_v<Range>);
		static_assert(parlay::is_less_than_comparable_v<
			      parlay::range_reference_type_t<Range>>);
		static_assert(std::is_constructible_v<
			      parlay::range_value_type_t<Range>,
			      parlay::range_reference_type_t<Range>>);

		auto aux = points_type::uninitialized(parlay::size(in));
		parlay::copy(in, parlay::make_slice(aux));
		op(parlay::make_slice(aux));
	}

	// NOTE: build tree
	static inline void sample_points(slice_type in, points_type &arr);

	static inline bucket_type
	find_bucket(Point const &p, hyper_plane_seq_type const &pivots);

	template <is_binary_node interior_type, bool UpdateParFlag = true>
	static inline void update_interior(node *T, node *L, node *R);

	template <is_binary_node interior_type, bool UpdateParFlag = true>
	static inline void update_interior(node *T, node_box_type const &L,
					   node_box_type const &R);

	template <is_multi_node interior_type, bool UpdateParFlag = true>
	static inline void
	update_interior(node *T,
			typename interior_type::node_arr_type const &new_nodes);

	template <typename leaf_type, typename interior_type,
		  bool granularity = true>
	static void prepare_rebuild(node *T, slice_type in, points_type &wx,
				    points_type &wo);

	template <typename leaf_type, typename interior_type,
		  bool granularity = true>
	static void prepare_rebuild(node *T, points_type &wx, points_type &wo);

	static void partition(slice_type A, slice_type B, size_t const n,
			      hyper_plane_seq_type const &pivots,
			      parlay::sequence<balls_type> &sums);

	// NOTE: batch insert
	template <typename leaf_type>
	static node *insert_points2_leaf(node *T, slice_type in);

	template <typename leaf_type, typename RT>
	static RT delete_points4_leaf(node *T, slice_type in);

	template <typename leaf_type, typename RT>
	static RT diff_points4_leaf(node *T, slice_type in);

	template <is_binary_node interior_type>
	static inline bucket_type retrieve_tag(Point const &p,
					       node_tag_seq_type const &tags);

	template <is_multi_node interior_type>
	static inline bucket_type retrieve_tag(Point const &p,
					       node_tag_seq_type const &tags);

	template <typename interior_type>
	static void sieve_points(slice_type A, slice_type B, size_t const n,
				 node_tag_seq_type const &tags,
				 parlay::sequence<balls_type> &sums,
				 bucket_type const tags_num);

	template <typename leaf_type, typename interior_type,
		  typename PrepareFunc, typename... Args>
	node *rebuild_with_insert(node *T, PrepareFunc prepare_func,
				  slice_type in, Args &&...args);

	template <typename leaf_type, typename interior_type, bool granularity,
		  typename... Args>
	node *rebuild_single_tree(node *T, Args &&...args);

	template <typename leaf_type, is_binary_node interior_type,
		  bool granularity, typename PrepareFunc, typename... Args>
	node *rebuild_tree_recursive(node *T, PrepareFunc &&prepare_func,
				     bool const allow_rebuild, Args &&...args);

	template <typename leaf_type, is_multi_node interior_type,
		  bool granularity, typename PrepareFunc, typename... Args>
	node *rebuild_tree_recursive(node *T, PrepareFunc &&prepare_func,
				     bool const allow_rebuild, Args &&...args);

	template <typename leaf_type, is_binary_node interior_type,
		  typename Base>
	static Base build_inner_tree(bucket_type idx,
				     hyper_plane_seq_type const &pivots,
				     parlay::sequence<Base> const &tree_nodes);

	template <is_multi_node interior_type>
	static node *
	build_inner_tree(bucket_type idx, hyper_plane_seq_type const &pivots,
			 parlay::sequence<node *> const &tree_nodes);

	// NOTE: delete tree
	template <supports_force_parallel interior_type, bool granularity>
	inline static bool force_parallel_recursion(interior_type const *T);

	/*
	 * Needed because delete_tree is virtual. The derived destructor does
	 * the freeing; calling a pure virtual from here would be too late.
	 */
	virtual ~base_tree() = default;

	base_tree() = default;

	/*
	 * root_ is an owning raw pointer, so a copy would be freed twice.
	 * A tree moves instead; the source is left empty and safe to destroy.
	 */
	base_tree(base_tree const &) = delete;
	base_tree &operator=(base_tree const &) = delete;

	base_tree(base_tree &&other) noexcept
	    : root_(std::exchange(other.root_, nullptr)),
	      tree_box_(std::exchange(other.tree_box_, get_empty_box()))
	{
	}

	base_tree &operator=(base_tree &&other) noexcept
	{
		if (this != &other) {
			delete_tree();
			root_ = std::exchange(other.root_, nullptr);
			tree_box_ =
				std::exchange(other.tree_box_, get_empty_box());
		}
		return *this;
	}

	constexpr virtual void delete_tree() = 0;

	template <typename leaf_type, typename interior_type>
	void delete_tree_wrapper();

	/*
	 * Frees the nodes and leaves the box alone, which is what build_ needs:
	 * delete_tree() would also drop a box the caller pinned.
	 */
	template <typename leaf_type, typename interior_type>
	void delete_tree_nodes();

	template <typename leaf_type, is_binary_node interior_type,
		  bool granularity = true>
	static void delete_tree_recursive(node *T);

	template <typename leaf_type, is_multi_node interior_type,
		  bool granularity = true>
	static void delete_tree_recursive(node *T);

	template <typename leaf_type, is_dynamic_node interior_type,
		  bool granularity = true>
	static void delete_tree_recursive(node *T);

	// NOTE: knn query stuffs
	static inline dis_type p2p_distance_square(Point const &p,
						   Point const &q);

	static inline dis_type p2b_min_distance_square(Point const &p,
						       box_type const &a);

	static inline dis_type p2b_max_distance_square(Point const &p,
						       box_type const &a);

	static inline double
	p2c_min_distance(Point const &p, Point const &center, dis_type const r);

	template <typename CircleType>
	static inline double p2c_min_distance(Point const &p,
					      CircleType const &cl);

	static inline dis_type
	interruptible_distance(Point const &p, Point const &q, dis_type up);

	// NOTE: searech knn in the leaf
	template <typename leaf_type, typename Range>
	static void knn_leaf(node *T, Point const &q,
			     bounded_queue<Point, Range> &bq);

	// NOTE: search knn in the binary node
	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void knn_binary(node *T, Point const &q,
			       bounded_queue<Point, Range> &bq,
			       box_type const &node_box, knn_logger &logger)
		requires std::same_as<typename interior_type::st_type,
				      hyper_plane_type>;

	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void knn_binary_box(node *T, Point const &q,
				   bounded_queue<Point, Range> &bq,
				   knn_logger &logger);

	// NOTE: search knn in the expanded multi node
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void
	knn_multi_expand(node *T, Point const &q, dims_type dim,
			 bucket_type idx, bounded_queue<Point, Range> &bq,
			 box_type const &node_box, knn_logger &logger);

	// NOTE: search knn in the expanded multi node
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void knn_multi_expand_box(node *T, Point const &q, dims_type dim,
					 bucket_type idx,
					 bounded_queue<Point, Range> &bq,
					 knn_logger &logger);

	// NOTE: search knn in the multi node
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void knn_multi(node *T, Point const &q,
			      bounded_queue<Point, Range> &bq,
			      knn_logger &logger);

	// NOTE: search knn in the cover node
	// TODO: change type of interior
	template <typename leaf_type, is_dynamic_node interior_type,
		  typename Range>
	static void knn_cover(node *T, Point const &q,
			      bounded_queue<Point, Range> &bq,
			      knn_logger &logger);

	// NOTE: range count stuffs
	template <typename leaf_type>
	static size_t range_count_rectangle_leaf(node *T,
						 box_type const &query_box);

	template <typename leaf_type, is_binary_node interior_type>
	static size_t range_count_rectangle(node *T, box_type const &query_box,
					    box_type const &node_box,
					    range_query_logger &logger);

	template <typename leaf_type, is_multi_node interior_type>
	static size_t range_count_rectangle(node *T, box_type const &query_box,
					    box_type const &node_box,
					    dims_type dim, bucket_type idx,
					    range_query_logger &logger);

	// template <typename leaf_type, is_binary_node interior_type>
	// static size_t RangeCountRadius(node* T, normal_circle const& cl,
	//                                box_type const& node_box);

	// NOTE: range query stuffs
	template <typename leaf_type, typename Range>
	static void range_query_leaf(node *T, Range out, size_t &s,
				     box_type const &query_box);

	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void range_query_serial_recursive(node *T, Range out, size_t &s,
						 box_type const &query_box,
						 box_type const &node_box,
						 range_query_logger &logger);

	// template <typename leaf_type, is_multi_node interior_type>
	// static size_t RangeCountRadius(node* T, normal_circle const& cl,
	//                                box_type const& node_box);

	// NOTE: range query stuffs
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void range_query_serial_recursive(node *T, Range out, size_t &s,
						 box_type const &query_box,
						 box_type const &node_box,
						 dims_type dim, bucket_type idx,
						 range_query_logger &logger);

	// NOTE: utility
	// TODO: better evaluate the parallel recursion function
	template <typename leaf_type, is_binary_node interior_type,
		  typename Range, bool granularity = true>
	static void flatten_rec(node *T, Range out);

	template <typename leaf_type, typename interior_type, typename Range,
		  bool granularity = true>
	static void flatten_rec(node *T, Range out)
		requires(!is_binary_node<interior_type>);

	template <typename leaf_type, typename Range>
	static void extract_points_in_leaf(node *T, Range out);

	template <typename leaf_type, is_multi_node interior_type,
		  typename Range, bool granularity = true>
	static void partial_flatten(node *T, Range out, bucket_type idx);

	// NOTE: validations
	template <typename leaf_type, typename interior_type>
	box_type check_box(node *T, box_type const &box);

	template <typename leaf_type, typename interior_type>
	points_type check_cover(
		node *T,
		typename interior_type::CircleType const &level_cover_circle);

	template <typename leaf_type, typename interior_type>
	static size_t check_size(node *T);

	template <typename leaf_type, typename interior_type>
	void check_tree_same_sequential(node *T, int dim);

	template <typename leaf_type, typename interior_type,
		  typename SplitRule>
	void validate();

	template <typename leaf_type, typename interior_type>
	size_t get_tree_height();

	template <typename leaf_type, typename interior_type>
	size_t get_max_tree_depth(node *T, size_t deep);

	template <typename leaf_type, typename interior_type>
	double get_ave_tree_height();

	template <typename leaf_type, typename interior_type>
	size_t count_tree_nodes_num(node *T);

	template <typename leaf_type, typename interior_type>
	void count_tree_heights(node *T, size_t deep, size_t &leaf_num,
				size_t &depth_sum);

	// NOTE: param interfaces
	void set_root(node *root)
	{
		this->root_ = root;
	}

	node *get_root() const
	{
		return this->root_;
	}

	size_t get_size() const
	{
		return this->root_ ? this->root_->size : 0;
	}

	box_type get_root_box() const
	{
		return this->tree_box_;
	}

	consteval static auto get_build_depth_once()
	{
		return static_cast<int>(build_depth_once);
	}

protected:
	node *root_ = nullptr;
	parlay::internal::timer timer;
	/* Empty box is the identity for get_box(); basic_point() leaves
	 * junk. */
	box_type tree_box_ = get_empty_box();
};

} // namespace psi

#include "psi/base_tree_impl/box_cut.hpp"
#include "psi/base_tree_impl/box_op.hpp"
#include "psi/base_tree_impl/circle_op.hpp"
#include "psi/base_tree_impl/delete_tree.hpp"
#include "psi/base_tree_impl/dimensionality.hpp"
#include "psi/base_tree_impl/inner_tree.hpp"
#include "psi/base_tree_impl/knn_query.hpp"
#include "psi/base_tree_impl/points_op.hpp"
#include "psi/base_tree_impl/range_query.hpp"
#include "psi/base_tree_impl/tree_op/build_inner_tree.hpp"
#include "psi/base_tree_impl/tree_op/flatten.hpp"
#include "psi/base_tree_impl/tree_op/leaf_op.hpp"
#include "psi/base_tree_impl/tree_op/node_op.hpp"
#include "psi/base_tree_impl/tree_op/rebuild.hpp"
#include "psi/base_tree_impl/validation.hpp"

#endif // PSI_BASE_TREE_H_
