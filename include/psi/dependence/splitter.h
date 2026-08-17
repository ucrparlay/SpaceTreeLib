#ifndef PSI_DEPENDENCE_SPLITTER_H
#define PSI_DEPENDENCE_SPLITTER_H

#include <concepts>
#include <cstdint>

#include "psi/tree_traits.h"
#include "psi/dependence/concepts.h"
#include "psi/dependence/tree_node.h"

#include "libmorton/morton.h"
#include "libmorton/morton3D.h"
#include "psi/dependence/space_filling_curve/hilbert.h"
#include "psi/dependence/space_filling_curve/hilbert_high_dim.h"

/*
 * hilbert.hpp holds the definitions, not just declarations, and is third
 * party. Its functions are inline so two translation units can both include
 * this header; its warnings are not ours to fix.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "psi/dependence/space_filling_curve/hilbert.hpp"
#pragma GCC diagnostic pop

namespace psi
{
/* ---------------- Spatial filling curve --------------- */
template <typename Point>
struct morton_curve {
	using base_type = point_traits<Point>;
	using slice_type = base_type::slice_type;
	using dims_type = base_type::dims_type;
	using box_type = base_type::box_type;
	using num_type = base_type::num_type;
	using coord_type = base_type::coord_type;
	using hyper_plane_type = base_type::hyper_plane_type;

	using curve_code_type = typename Point::at_type::curve_code_type;

	void morton_tag()
	{
	}
	static std::string get_name()
	{
		return "MortonCurve";
	}

	static auto encode(Point const &p)
	{
		assert(std::is_integral_v<coord_type>);
		if constexpr (Point::get_dim() == 2) {
			return libmorton::m2D_e_for<curve_code_type, uint32_t>(
				p.pnt[0], p.pnt[1]);
		} else if constexpr (Point::get_dim() == 3) {
			return libmorton::m3D_e_for_ET<curve_code_type,
						       uint32_t>(
				p.pnt[0], p.pnt[1], p.pnt[2]);
		} else {
			static_assert(
				"MortonCurve only supports 2D and 3D points");
		}
	}
};

template <typename Point>
struct hilbert_curve {
	using base_type = point_traits<Point>;
	using slice_type = base_type::slice_type;
	using dims_type = base_type::dims_type;
	using box_type = base_type::box_type;
	using num_type = base_type::num_type;
	using coord_type = base_type::coord_type;
	using hyper_plane_type = base_type::hyper_plane_type;

	using curve_code_type = typename Point::at_type::curve_code_type;

	void hilbert_tag()
	{
	}
	static std::string get_name()
	{
		return "HilbertCurve";
	}

	/* TODO: optimize the encode */
	static auto encode(Point const &p)
	{
		assert(std::is_integral_v<coord_type>);

		if constexpr (Point::get_dim() == 2) {
			return hilbert::hilbert_c2i(
				2, 32,
				reinterpret_cast<hilbert::bitmask_t const *>(
					p.get_coords().data()));
		} else if constexpr (Point::get_dim() == 3) {
			/* maximum 20 bits for each dimension (2^20 = 1048576))
			 */
			return hilbert::hilbert_c2i(
				3, 21,
				reinterpret_cast<hilbert::bitmask_t const *>(
					p.get_coords().data()));
		} else {
			static_assert(
				"HilbertCurve only supports 2D and 3D points");
		}
	}
};

template <typename Curve>
struct spatial_filling_curve {
	using curve_code_type = typename Curve::curve_code_type;

	static std::string get_split_name()
	{
		return Curve::get_name();
	}

	template <typename... Args>
	static auto encode(Args &&...args)
	{
		return Curve::encode(std::forward<Args>(args)...);
	}

	template <typename... Args>
	auto decode(Args &&...args)
	{
		return Curve::decode(std::forward<Args>(args)...);
	}
};

/* ---------------- Orthogonal Split Rule ---------------- */
template <typename Point>
struct base_split_dim_rule {
	using base_type = point_traits<Point>;
	using slice_type = base_type::slice_type;
	using dims_type = base_type::dims_type;
	using box_type = base_type::box_type;
	using num_type = base_type::num_type;
	using coord_type = base_type::coord_type;
	using hyper_plane_type = base_type::hyper_plane_type;

	constexpr virtual dims_type const
	find_cutting_dimension(box_type const &bx,
			       dims_type const dim) const = 0;

	constexpr virtual dims_type const
	find_rebuild_dimension(dims_type const dim) const = 0;

	constexpr virtual dims_type const
	next_dimension(dims_type const dim) const = 0;
	/* TODO: the spliiter should deterine how to split as well */
};

template <typename Point>
struct max_stretch_dim : base_split_dim_rule<Point> {
	using bsr_type = base_split_dim_rule<Point>;
	using base_type = bsr_type::base_type;
	using slice_type = bsr_type::slice_type;
	using dims_type = bsr_type::dims_type;
	using box_type = bsr_type::box_type;
	using coord_type = bsr_type::coord_type;
	using num_type = bsr_type::num_type;
	using hyper_plane_type = base_type::hyper_plane_type;

	void max_stretch_tag()
	{
	}
	static std::string get_name()
	{
		return "MaxStretchDim";
	}

	constexpr dims_type const
	find_cutting_dimension(box_type const &bx, dims_type) const override
	{
		dims_type d(0);
		coord_type diff(bx.second.pnt[0] - bx.first.pnt[0]);
		assert(num_type::geq(diff, 0));
		for (dims_type i = 1; i < base_type::num_dims; ++i) {
			if (num_type::gt(bx.second.pnt[i] - bx.first.pnt[i],
					 diff)) {
				diff = bx.second.pnt[i] - bx.first.pnt[i];
				d = i;
			}
		}
		return d;
	};

	constexpr dims_type const
	find_rebuild_dimension(dims_type) const override
	{
		return 0;
	};

	/*
	 * TODO: this is wired as the next dimension should return the dimension
	 * with second largest
	 */
	constexpr dims_type const next_dimension(dims_type) const override
	{
		return 0;
	};
};

template <typename Point>
struct rotate_dim : base_split_dim_rule<Point> {
	using bsr_type = base_split_dim_rule<Point>;
	using base_type = bsr_type::base_type;
	using slice_type = bsr_type::slice_type;
	using dims_type = bsr_type::dims_type;
	using box_type = bsr_type::box_type;
	using coord_type = bsr_type::coord_type;
	using num_type = bsr_type::num_type;
	using points_iter_type = base_type::points_iter_type;
	using hyper_plane_type = base_type::hyper_plane_type;

	void rotate_dim_tag()
	{
	}
	static std::string get_name()
	{
		return "RotateDim";
	}

	constexpr dims_type const
	find_cutting_dimension([[maybe_unused]] box_type const &bx,
			       dims_type const dim) const override
	{
		return dim;
	};

	constexpr dims_type const find_rebuild_dimension(
		[[maybe_unused]] dims_type const dim) const override
	{
		return dim;
	};

	constexpr dims_type const next_dimension(dims_type dim) const override
	{
		return (dim + 1) % base_type::num_dims;
	};
};

template <typename Point>
struct base_split_partition_rule {
	using base_type = point_traits<Point>;
	using slice_type = base_type::slice_type;
	using dims_type = base_type::dims_type;
	using box_type = base_type::box_type;
	using num_type = base_type::num_type;
	using coord_type = base_type::coord_type;
	using hyper_plane_type = base_type::hyper_plane_type;
	using points_iter_type = base_type::points_iter_type;
	using iter_hyper_pair_type =
		std::pair<points_iter_type, std::optional<hyper_plane_type>>;

	constexpr virtual iter_hyper_pair_type const
	split_input(slice_type in, dims_type const dim,
		    box_type const &box) const = 0;
	constexpr virtual hyper_plane_type const
	split_sample(slice_type in, dims_type const dim,
		     box_type const &box) const = 0;

	/*
	 * Batch update asks this before rebuilding a subtree. It was defined
	 * ad hoc on each concrete rule and absent here, so a new rule that
	 * omitted it built cleanly and failed on the first batch update
	 * instead of at the declaration.
	 */
	constexpr virtual bool allow_rebuild() const = 0;
};

template <typename Point>
struct object_median : base_split_partition_rule<Point> {
	using bsr_type = base_split_partition_rule<Point>;
	using base_type = bsr_type::base_type;
	using slice_type = bsr_type::slice_type;
	using dims_type = bsr_type::dims_type;
	using box_type = bsr_type::box_type;
	using coord_type = bsr_type::coord_type;
	using num_type = bsr_type::num_type;
	using points_iter_type = bsr_type::points_iter_type;
	using hyper_plane_type = bsr_type::hyper_plane_type;
	using iter_hyper_pair_type = bsr_type::iter_hyper_pair_type;

	void object_median_tag()
	{
	}

	static std::string get_name()
	{
		return "ObjectMedian";
	}

	constexpr bool allow_rebuild() const override
	{
		return true;
	};

	constexpr iter_hyper_pair_type const
	split_input(slice_type in, dims_type const dim,
		    [[maybe_unused]] box_type const &box) const override
	{
		size_t n = in.size();
		std::ranges::nth_element(
			in.begin(), in.begin() + n / 2, in.end(),
			[&](Point const &p1, Point const &p2) {
				return num_type::lt(p1.pnt[dim], p2.pnt[dim]);
			});

		auto split_iter =
			std::ranges::partition(
				in.begin(), in.begin() + n / 2,
				[&](Point const &p) {
					return num_type::lt(p.pnt[dim],
							    in[n / 2].pnt[dim]);
				})
				.begin();

		if (split_iter == in.begin()) { /* handle duplicated medians */
			split_iter =
				std::ranges::partition(
					in.begin() + n / 2, in.end(),
					[&](Point const &p) {
						return num_type::eq(
							p.pnt[dim],
							in[n / 2].pnt[dim]);
					})
					.begin(); /* now all duplicated */
						  /* median is on the left */
		}
		if (split_iter <=
		    in.begin() + n / 2) { /* split is on left half */
			return iter_hyper_pair_type(
				split_iter,
				hyper_plane_type(in[n / 2].pnt[dim], dim));
		} else if (split_iter !=
			   in.end()) { /* split is on right half */
			auto min_elem_iter = std::ranges::min_element(
				split_iter, in.end(),
				[&](Point const &p1, Point const &p2) {
					return num_type::lt(p1.pnt[dim],
							    p2.pnt[dim]);
				});
			return iter_hyper_pair_type(
				split_iter,
				hyper_plane_type(min_elem_iter->pnt[dim], dim));
		} else { /* all the same */
			return iter_hyper_pair_type(split_iter, std::nullopt);
		}
	}

	/* TODO: the sample may handle the duplicates as well */
	constexpr hyper_plane_type const
	split_sample(slice_type in, dims_type const dim,
		     [[maybe_unused]] box_type const &box) const override
	{
		size_t n = in.size();
		std::ranges::nth_element(
			in, in.begin() + n / 2,
			[&](Point const &p1, Point const &p2) {
				return num_type::lt(p1.pnt[dim], p2.pnt[dim]);
			});
		return hyper_plane_type(in[n / 2][dim], dim);
	}
};

template <typename Point>
struct spatial_median : base_split_partition_rule<Point> {
	using bsr_type = base_split_partition_rule<Point>;
	using base_type = bsr_type::base_type;
	using slice_type = bsr_type::slice_type;
	using dims_type = bsr_type::dims_type;
	using box_type = bsr_type::box_type;
	using coord_type = bsr_type::coord_type;
	using num_type = bsr_type::num_type;
	using points_iter_type = bsr_type::points_iter_type;
	using hyper_plane_type = bsr_type::hyper_plane_type;
	using iter_hyper_pair_type = bsr_type::iter_hyper_pair_type;

	void spatial_median_tag()
	{
	}

	static std::string get_name()
	{
		return "SpatialMedian";
	}

	constexpr bool allow_rebuild() const override
	{
		return false;
	};

	constexpr iter_hyper_pair_type const
	split_input(slice_type in, dims_type const dim,
		    box_type const &box) const override
	{
		auto split_iter =
			std::ranges::partition(in, [&](Point const &p) {
				return num_type::lt(
					p.pnt[dim],
					base_type::get_box_mid(dim, box));
			}).begin();
		if (split_iter == in.begin() || split_iter == in.end()) {
			return iter_hyper_pair_type(split_iter, std::nullopt);
		} else {
			return iter_hyper_pair_type(
				split_iter,
				hyper_plane_type(
					base_type::get_box_mid(dim, box), dim));
		}
	}

	constexpr hyper_plane_type const
	split_sample([[maybe_unused]] slice_type in, dims_type const dim,
		     box_type const &box) const override
	{
		return hyper_plane_type(base_type::get_box_mid(dim, box), dim);
	}
};

template <class DimRule, class PartitionRule>
struct orthogonal_split_rule {
	using dim_rule_type = DimRule;
	using partition_rule_type = PartitionRule;

	static std::string get_split_name()
	{
		return DimRule::get_name() + "-" + PartitionRule::get_name();
	}

	/* dimension */
	template <typename... Args>
	auto find_cutting_dimension(Args &&...args)
	{
		return dim_rule.find_cutting_dimension(
			std::forward<Args>(args)...);
	}

	template <typename... Args>
	auto find_rebuild_dimension(Args &&...args)
	{
		return dim_rule.find_rebuild_dimension(
			std::forward<Args>(args)...);
	}

	template <typename... Args>
	auto next_dimension(Args &&...args)
	{
		return dim_rule.next_dimension(std::forward<Args>(args)...);
	}

	/* serial parititon used in algorithm */
	template <typename... Args>
	auto split_input(Args &&...args)
	{
		return partition_rule.split_input(std::forward<Args>(args)...);
	}

	/* split the sample in order to get the hyperplane */
	template <typename... Args>
	auto split_sample(Args &&...args)
	{
		return partition_rule.split_sample(std::forward<Args>(args)...);
	}

	/* query whether to launch the rebuild */
	template <typename... Args>
	auto allow_rebuild(Args &&...args)
	{
		return partition_rule.allow_rebuild(
			std::forward<Args>(args)...);
	}

	/*
	 * helper for handling the duplicate
	 * divide the space until the split cut the input box
	 * INFO: divide the space for binary node
	 */
	template <typename Tree, typename slice_type, typename dims_type,
		  typename box_type>
	node *divide_space(Tree &tree, slice_type in, slice_type out,
			   box_type const &node_box, box_type const &input_box,
			   dims_type dim)
		requires(is_binary_node<typename Tree::interior_type>)
	{
		assert(Tree::within_box(input_box, node_box));

		/* Main logic */
		if (Tree::vertical_line_split_box(
			    Tree::get_box_mid(dim, node_box), input_box, dim)) {
			return tree.build_recursive(in, out, dim, node_box);
		}

		auto cut_dim = dim_rule.find_cutting_dimension(node_box, dim);
		auto cut_val = Tree::get_box_mid(cut_dim, node_box);

		/*
		 * INFO: Test whether the node_box will remain unchanged after
		 * the split. If the mid of box is the same as the box edge,
		 * then this time the recursion will usless, the worst case is
		 * that all the mid on all dimension is on the box edge, i.e.,
		 * (0,1), (0,1), then a correct split algorithm will handle this
		 * case
		 */
		dims_type dim_cnt = 0;
		while (dim_cnt != Tree::num_dims) {
			if (!Tree::vertical_line_on_box_edge(cut_val, node_box,
							     cut_dim)) {
				break;
			}
			cut_dim = dim_rule.next_dimension(cut_dim);
			cut_val = Tree::get_box_mid(cut_dim, node_box);
			dim_cnt++;
		}

		if (dim_cnt ==
		    Tree::num_dims) { /* this breaks rotation manner */
			/* checks whether the node box is separatable */
			return tree.build_recursive(
				in, out, dim_rule.next_dimension(dim),
				node_box);
		} else if (Tree::vertical_line_split_box(cut_val, input_box,
							 cut_dim)) {
			/*
			 * above while loop may changed to new dim and
			 * need to re-check whether it can split the input. This
			 * is necessary as the following left/right test assumes
			 * the @cut_val does not split box
			 */
			return tree.build_recursive(in, out, cut_dim, node_box);
		}

		bool split_is_right =
			Tree::num_type::gt(cut_val, input_box.second[cut_dim]);
		assert(split_is_right ||
		       Tree::num_type::leq(cut_val, input_box.first[cut_dim]));

		typename Tree::box_cut_type box_cut(
			node_box,
			typename Tree::splitter_type(cut_val, cut_dim),
			split_is_right);

		node *L = divide_space(tree, in, out,
				       box_cut.get_first_box_cut(), input_box,
				       dim_rule.next_dimension(cut_dim));
		node *R = alloc_empty_leaf_node<slice_type,
						typename Tree::leaf_type>();
		assert(Tree::within_box(input_box, box_cut.get_box()));

		if (!split_is_right) {
			assert(Tree::num_type::leq(
				Tree::get_box_mid(cut_dim, node_box),
				input_box.first[cut_dim]));
			std::ranges::swap(L, R);
		}

		return alloc_interior_node<typename Tree::interior_type>(
			L, R, box_cut.get_hyper_plane());
	}

	/* INFO: divide the space for multi node */
	template <typename Tree, typename slice_type, typename box_type>
	node *divide_space(Tree &tree, slice_type in, slice_type out,
			   box_type const &node_box, box_type const &input_box)
		requires(is_multi_node<typename Tree::interior_type>)
	{
		assert(Tree::within_box(input_box, node_box));

		auto nodebox_split =
			Tree::interior_type::compute_splitter(node_box);
		for (auto const &split : nodebox_split) {
			if (Tree::vertical_line_split_box(
				    split.first, input_box, split.second)) {
				/*
				 * any point on the split should be put on the
				 * right
				 */
				return tree.build_recursive(in, out, node_box);
			}
		}

		if (std::ranges::all_of(nodebox_split, [&](auto const &split) {
			    return Tree::vertical_line_on_box_edge(
				    split.first, node_box, split.second);
		    })) { /* new box will be same as the old ones */
			return tree.build_recursive(in, out, node_box);
		}

		/* PARA: node id contains the input box */
		typename Tree::interior_type::bucket_type pivot_region_id = 1;
		/*
		 * considering three cases:
		 * 1: the split is LESS than the left input box edge
		 * 2: the split is EQ the left box edge
		 * 3: same as 2 but left and right edge is same, aka. a line
		 * The special judge above has prune out the case that split
		 * intersect the box, therefore split lt or eq the RIGHT edge
		 * will cover all the cases
		 */
		for (auto const &split : nodebox_split) {
			pivot_region_id =
				pivot_region_id * 2 +
				static_cast<typename Tree::interior_type::
						    bucket_type>(
					Tree::num_type::leq(
						split.first,
						input_box
							.second[split.second]));
		}
		pivot_region_id -= Tree::interior_type::get_regions();

		typename Tree::interior_type::node_arr_type tree_nodes;
		for (typename Tree::interior_type::bucket_type i = 0;
		     i < Tree::interior_type::get_regions(); i++) {
			tree_nodes[i] =
				pivot_region_id == i
					? divide_space(
						  tree, in, out,
						  Tree::interior_type::
							  get_box_by_region_id(
								  i,
								  nodebox_split,
								  node_box),
						  input_box)
					: alloc_empty_leaf_node<
						  slice_type,
						  typename Tree::leaf_type>();
		}

		return alloc_interior_node<typename Tree::interior_type>(
			tree_nodes, nodebox_split);
	}

	/*
	 * cannot divide the points on @dim, while the points are not the
	 * same
	 */
	template <typename Tree, typename slice_type, typename box_type,
		  typename... Args>
	auto handle_undivided(Tree &tree, slice_type in, slice_type out,
			      box_type const &box, Args &&...args)
	{
		if constexpr (is_object_median_split<PartitionRule>) {
			/*
			 * in object median, if current dimension is not
			 * divideable, then switch to another dimension then
			 * continue. This works since unless all points are
			 * same, otherwise we can always slice some points out.
			 */
			if constexpr (is_binary_node<
					      typename Tree::interior_type>) {
				/* TODO: tooo brute force */
				return tree.serial_build_recursive(
					in, out,
					dim_rule.next_dimension(
						std::forward<Args>(args)...),
					Tree::get_box(in));
			} else {
				return tree.build_recursive(in, out,
							    Tree::get_box(in));
			}

		} else if constexpr (is_spatial_median_split<PartitionRule>) {
			/*
			 * in spatial median, we simply reduce the box by
			 * half on current dim, then switch to next dim.
			 */
			if constexpr (is_rotate_dim_split<DimRule> ||
				      is_max_stretch_dim<DimRule>) {
				return divide_space(
					tree, in, out, box, Tree::get_box(in),
					std::forward<Args>(args)...);
			} else {
				/*
				 * A rule with no branch here falls out of a
				 * function that has to return a node, and the
				 * contributor sees "void value not ignored"
				 * in a file they never opened. Say it where
				 * the edit belongs. Dependent, so it only
				 * fires on instantiation.
				 */
				static_assert(sizeof(DimRule) == 0,
					      "handle_undivided has no branch "
					      "for this dimension rule");
			}
		} else {
			static_assert(sizeof(PartitionRule) == 0,
				      "handle_undivided has no branch for "
				      "this partition rule");
		}
	}

	DimRule const dim_rule;
	PartitionRule const partition_rule;
};
} /* namespace psi */

#endif /* PSI_DEPENDENCE_SPLITTER_H */
