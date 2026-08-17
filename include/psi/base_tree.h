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
#include "psi/tree_traits.h"

namespace psi
{

/* Basetree */
template <typename Traits, typename DerivedTree>
class base_tree : public Traits
{
public:
	using traits_type = Traits;

	/*
	 * Traits is a dependent base, so its names are not visible unqualified
	 * inside the tree or in any of the impl files. Everything a tree needs
	 * is pulled in once, here.
	 */
	using point_type = typename Traits::point_type;
	using basic_point = typename Traits::basic_point;
	using coord_type = typename Traits::coord_type;
	using coords_type = typename Traits::coords_type;
	using dis_type = typename Traits::dis_type;
	using num_type = typename Traits::num_type;

	using bucket_type = typename Traits::bucket_type;
	using balls_type = typename Traits::balls_type;
	using dims_type = typename Traits::dims_type;
	using depth_type = typename Traits::depth_type;
	using bucket_seq_type = typename Traits::bucket_seq_type;
	using ball_seq_type = typename Traits::ball_seq_type;

	using slice_type = typename Traits::slice_type;
	using points_type = typename Traits::points_type;
	using points_iter_type = typename Traits::points_iter_type;
	using hyper_plane_type = typename Traits::hyper_plane_type;
	using hyper_plane_seq_type = typename Traits::hyper_plane_seq_type;
	using box_type = typename Traits::box_type;
	using box_seq_type = typename Traits::box_seq_type;
	using knn_result_type = typename Traits::knn_result_type;
	using knn_result_seq_type = typename Traits::knn_result_seq_type;

	using node_box_type = typename Traits::node_box_type;
	using node_box_seq_type = typename Traits::node_box_seq_type;
	using node_tag_type = typename Traits::node_tag_type;
	using node_tag_seq_type = typename Traits::node_tag_seq_type;

	using split_rule_type = typename Traits::split_rule_type;

	using Traits::block_size;
	using Traits::bucket_num;
	using Traits::build_depth_once;
	using Traits::imbalance_ratio;
	using Traits::leaf_capacity;
	using Traits::log2_base;
	using Traits::num_dims;
	using Traits::pivot_num;
	using Traits::serial_build_cutoff;
	using Traits::sparse_leaf_threshold;

	using Traits::box_intersect_box;
	using Traits::box_is_empty;
	using Traits::get_box;
	using Traits::get_box_center;
	using Traits::get_box_mid;
	using Traits::get_empty_box;
	using Traits::is_box_line_in_dimension;
	using Traits::legal_box;
	using Traits::same_box;
	using Traits::vertical_line_intersect_box;
	using Traits::vertical_line_intersect_box_exclude;
	using Traits::vertical_line_on_box_edge;
	using Traits::vertical_line_on_box_left_edge;
	using Traits::vertical_line_on_box_right_edge;
	using Traits::vertical_line_split_box;
	using Traits::within_box;

	/* array based inner tree for batch insertion and deletion */
	template <typename leaf_type, typename interior_type>
	struct inner_tree;

	/* compute the bounding box on the fly */
	struct box_cut_type;

	/* get the imbalance ratio */
	static inline size_t get_imbalance_ratio();
	static inline bool imbalance_node(size_t const l, size_t const n);
	static inline bool sparse_node(size_t const l, size_t const n);

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

	/* build tree */
	static inline void sample_points(slice_type in, points_type &arr);

	static inline bucket_type
	find_bucket(point_type const &p, hyper_plane_seq_type const &pivots);

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

	/* batch insert */
	template <typename leaf_type>
	static node *insert_points2_leaf(node *T, slice_type in);

	template <typename leaf_type, typename RT>
	static RT delete_points4_leaf(node *T, slice_type in);

	template <typename leaf_type, typename RT>
	static RT diff_points4_leaf(node *T, slice_type in);

	template <is_binary_node interior_type>
	static inline bucket_type retrieve_tag(point_type const &p,
					       node_tag_seq_type const &tags);

	template <is_multi_node interior_type>
	static inline bucket_type retrieve_tag(point_type const &p,
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

	/* delete tree */
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

	/*
	 * Convenience forms. The primary entry points write into a buffer the
	 * caller has to size, which means knowing the answer first; these
	 * allocate and return. They are ordinary wrappers -- nothing here is
	 * on a measured path unless a caller chooses it.
	 */
	points_type range_query(box_type const &query_box) const;
	knn_result_seq_type knn(point_type const &q, size_t k) const;

	/* Through the derived tree: p_tree keeps its points in a cpam map and
	 * has no root_ at all. */
	bool empty() const
	{
		return static_cast<DerivedTree const *>(this)->get_size() == 0;
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

	/* knn query stuffs */
	static inline dis_type p2p_distance_square(point_type const &p,
						   point_type const &q);

	static inline dis_type p2b_min_distance_square(point_type const &p,
						       box_type const &a);

	static inline dis_type p2b_max_distance_square(point_type const &p,
						       box_type const &a);

	static inline double p2c_min_distance(point_type const &p,
					      point_type const &center,
					      dis_type const r);

	template <typename CircleType>
	static inline double p2c_min_distance(point_type const &p,
					      CircleType const &cl);

	static inline dis_type interruptible_distance(point_type const &p,
						      point_type const &q,
						      dis_type up);

	/* searech knn in the leaf */
	template <typename leaf_type, typename Range>
	static void knn_leaf(node *T, point_type const &q,
			     bounded_queue<point_type, Range> &bq);

	/* search knn in the binary node */
	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void knn_binary(node *T, point_type const &q,
			       bounded_queue<point_type, Range> &bq,
			       box_type const &node_box, knn_logger &logger)
		requires std::same_as<typename interior_type::st_type,
				      hyper_plane_type>;

	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void knn_binary_box(node *T, point_type const &q,
				   bounded_queue<point_type, Range> &bq,
				   knn_logger &logger);

	/* search knn in the expanded multi node */
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void
	knn_multi_expand(node *T, point_type const &q, dims_type dim,
			 bucket_type idx, bounded_queue<point_type, Range> &bq,
			 box_type const &node_box, knn_logger &logger);

	/* search knn in the expanded multi node */
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void knn_multi_expand_box(node *T, point_type const &q,
					 dims_type dim, bucket_type idx,
					 bounded_queue<point_type, Range> &bq,
					 knn_logger &logger);

	/* search knn in the multi node */
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void knn_multi(node *T, point_type const &q,
			      bounded_queue<point_type, Range> &bq,
			      knn_logger &logger);

	/* range count stuffs */
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

	/* range query stuffs */
	template <typename leaf_type, typename Range>
	static void range_query_leaf(node *T, Range out, size_t &s,
				     box_type const &query_box);

	template <typename leaf_type, is_binary_node interior_type,
		  typename Range>
	static void range_query_serial_recursive(node *T, Range out, size_t &s,
						 box_type const &query_box,
						 box_type const &node_box,
						 range_query_logger &logger);

	/* range query stuffs */
	template <typename leaf_type, is_multi_node interior_type,
		  typename Range>
	static void range_query_serial_recursive(node *T, Range out, size_t &s,
						 box_type const &query_box,
						 box_type const &node_box,
						 dims_type dim, bucket_type idx,
						 range_query_logger &logger);

	/*
	 * utility
	 * TODO: better evaluate the parallel recursion function
	 */
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

	/* validations */
	template <typename leaf_type, typename interior_type>
	box_type check_box(node *T, box_type const &box);

	template <typename leaf_type, typename interior_type>
	static size_t check_size(node *T);

	template <typename leaf_type, typename interior_type>
	void check_tree_same_sequential(node *T, int dim);

	template <typename leaf_type, typename interior_type>
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

	node *get_root() const
	{
		return this->root_;
	}

	size_t get_size() const
	{
		return this->root_ ? this->root_->size : 0;
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

} /* namespace psi */

#include "psi/base_tree_impl/box_cut.hpp"
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

#endif /* PSI_BASE_TREE_H_ */
