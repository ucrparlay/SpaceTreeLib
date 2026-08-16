#ifndef PSI_BASE_TREE_IMPL_VALIDATION_HPP_
#define PSI_BASE_TREE_IMPL_VALIDATION_HPP_

#include <parlay/parallel.h>

#include <concepts>

#include "../base_tree.h"
#include "parlay/primitives.h"

namespace psi
{

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::check_box(
	node *T, box_type const &box)
{
	if (T->is_leaf) {
		auto const *tl = static_cast<leaf_type const *>(T);
		auto const node_box = get_box<leaf_type, interior_type>(T);
		if constexpr (has_box<typename leaf_type::at_type>) { // whether
								      // has bb
			assert(same_box(node_box, tl->get_box()));
		} else {
			assert(within_box(node_box, box));
		}
		return node_box;
	}
	interior_type *ti = static_cast<interior_type *>(T);
	assert(!ti->get_parallel_flag_ini_status()); // NOTE: ensure that
						     // uninitialized force
						     // parallelism
	if constexpr (is_binary_node<interior_type> &&
		      !has_box<typename interior_type::at_type>) { // NOTE: use
								   // hyperplane
		// as splitter
		box_type lbox(box), rbox(box);
		lbox.second.pnt[ti->split.second] = ti->split.first;
		rbox.first.pnt[ti->split.second] = ti->split.first;
		box_type const left_return_box =
			check_box<leaf_type, interior_type>(ti->left, lbox);
		box_type const right_return_box =
			check_box<leaf_type, interior_type>(ti->right, rbox);
		box_type const new_box =
			get_box(left_return_box, right_return_box);

		assert(within_box(left_return_box, lbox));
		assert(within_box(right_return_box, rbox));
		assert(within_box(new_box, box));
		return new_box;
	} else if constexpr (is_binary_node<interior_type> &&
			     has_box<typename interior_type::at_type>) { // kd
									 // with
									 // box
		// std::cout << " has box " << std::endl;
		auto left_box =
			retrieve_box<leaf_type, interior_type>(ti->left);
		auto right_box =
			retrieve_box<leaf_type, interior_type>(ti->right);
		box_type const left_return_box =
			check_box<leaf_type, interior_type>(ti->left, left_box);
		box_type const right_return_box =
			check_box<leaf_type, interior_type>(ti->right,
							    right_box);
		box_type const new_box =
			get_box(left_return_box, right_return_box);
		assert(same_box(
			left_return_box,
			retrieve_box<leaf_type, interior_type>(ti->left)));
		assert(same_box(
			right_return_box,
			retrieve_box<leaf_type, interior_type>(ti->right)));
		assert(same_box(new_box, ti->get_box()));
		return new_box;
	} else if (is_multi_node<interior_type> &&
		   !has_box<typename interior_type::at_type>) { // orth without
								// box
		box_seq_type new_box(
			ti->template compute_subregions<box_seq_type>(box));
		box_seq_type return_box_seq(new_box.size());
		assert(new_box.size() == ti->tree_nodes.size());
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			return_box_seq[i] = check_box<leaf_type, interior_type>(
				ti->tree_nodes[i], new_box[i]);
			assert(within_box(return_box_seq[i], new_box[i]));
		}
		auto return_box = get_box(return_box_seq);
		assert(within_box(return_box, box));
		return return_box;
	} else if (is_multi_node<interior_type> &&
		   has_box<typename interior_type::at_type>) { // orth with box
		box_seq_type new_box(
			ti->template compute_subregions<box_seq_type>(box));
		box_seq_type return_box_seq(new_box.size());
		assert(new_box.size() == ti->tree_nodes.size());
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			return_box_seq[i] = check_box<leaf_type, interior_type>(
				ti->tree_nodes[i], new_box[i]);
			assert(same_box(return_box_seq[i],
					retrieve_box<leaf_type, interior_type>(
						ti->tree_nodes[i])));
		}
		auto return_box = get_box(return_box_seq);
		assert(same_box(return_box, ti->get_box()));
		return return_box;
	} else {
		static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>,
			      "check_box supports only binary and multi-way "
			      "interior nodes");
		return get_empty_box();
	}
}

// NOTE: 1: whether the subtree is within the circle
// 2: separation of nodes per level
// 3: decrease level per level
// 4: nesting
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::points_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::check_cover(
	node *T, typename interior_type::CircleType const &level_cover_circle)
{
	using tree_circle_type = typename interior_type::CircleType;

	if (T->is_leaf) {
		auto tl = static_cast<leaf_type *>(T);
		for (auto const &p : tl->get_points()) {
			assert(within_circle(p, level_cover_circle));
		}
		return tl->get_points();
	}

	auto ti = static_cast<interior_type *>(T);
	assert(level_cover_circle == ti->get_cover_circle());

	// NOTE: nesting
	assert(std::ranges::count(ti->split, level_cover_circle.center) == 1);

	// NOTE: separation
	assert(std::all_of(
		ti->split.begin(), ti->split.end(), [&](auto const &center_i) {
			auto i = &center_i - &ti->split[0];
			auto tmp_circle = cover_circle{
				center_i, level_cover_circle.level - 1};
			return std::all_of(
				ti->split.begin() + i + 1, ti->split.end(),
				[&](auto const &center_j) {
					// PERF: should be gt as a point in a
					// boundary should falls within that
					// circle
					return num_type::gt(
						p2p_distance_square(center_i,
								    center_j),
						tmp_circle.get_radius_square());
				});
		}));

	// NOTE: covering
	// the cover circle of all points within the subtree should within the
	// cover circle of the node
	parlay::sequence<points_type> return_points_seq(ti->tree_nodes.size());
	for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
		auto next_circle = tree_circle_type{
			ti->split[i],
			static_cast<depth_type>(ti->get_cover_circle().level -
						1)};
		return_points_seq[i] = check_cover<leaf_type, interior_type>(
			ti->tree_nodes[i], next_circle);
		parlay::all_of(return_points_seq[i], [&](auto const &p) {
			return within_circle(p, level_cover_circle);
		});
		// for (auto const& p : return_points_seq[i]) {
		//   assert(within_circle(p, next_circle));
		// }
	}
	return parlay::flatten(return_points_seq);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
size_t base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::check_size(node *T)
{
	if (T->is_leaf) {
		// assert(static_cast<leaf_type*>(T)->is_dummy ||
		//        static_cast<leaf_type*>(T)->size <= leaf_capacity);
		return T->size;
	}
	if constexpr (is_binary_node<interior_type>) {
		interior_type *ti = static_cast<interior_type *>(T);
		assert(!ti->get_parallel_flag_ini_status());
		size_t l, r;
		parlay::par_do(
			[&l, &ti] {
				l = check_size<leaf_type, interior_type>(
					ti->left);
			},
			[&r, &ti] {
				r = check_size<leaf_type, interior_type>(
					ti->right);
			});
		assert(l + r == T->size);
		return T->size;
	} else {
		// assert(is_multi_node<interior_type>);
		interior_type *ti = static_cast<interior_type *>(T);
		assert(!ti->get_parallel_flag_ini_status());
		size_t sum = 0;
		for (bucket_type i = 0; i < ti->tree_nodes.size(); ++i) {
			sum += check_size<leaf_type, interior_type>(
				ti->tree_nodes[i]);
		}
		assert(sum == T->size);
		return T->size;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
void base_tree<Point, DerivedTree, SkHeight,
	       ImbaRatio>::check_tree_same_sequential(node *T, int dim)
{
	if (T->is_leaf) {
		// assert( PickRebuildDim( T, num_dims ) == dim );
		return;
	}
	if constexpr (is_binary_node<interior_type>) {
		interior_type *ti = static_cast<interior_type *>(T);
		if (ti->split.second != dim) {
			std::cout << int(ti->split.second) << " " << int(dim)
				  << " " << ti->size;
		}
		assert(ti->split.second == dim);
		// TODO: maybe need to add the split rule in the base tree?
		dim = (dim + 1) % num_dims;
		parlay::par_do_if(
			T->size > 1000,
			[&]() {
				check_tree_same_sequential<leaf_type,
							   interior_type>(
					ti->left, dim);
			},
			[&]() {
				check_tree_same_sequential<leaf_type,
							   interior_type>(
					ti->right, dim);
			});
	} else {
		assert(is_multi_node<interior_type>);
		interior_type *ti = static_cast<interior_type *>(T);
		assert(std::cmp_equal((1 << ti->split.size()),
				      ti->tree_nodes.size()));
		for (size_t i = 0; i < ti->split.size(); i++) {
			assert(ti->split[i].second == dim);
			dim += 1;
		}
		assert(dim == num_dims);
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			check_tree_same_sequential<leaf_type, interior_type>(
				ti->tree_nodes[i], 0);
		}
	}
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, typename SplitRule>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::validate()
{
	std::cout << ">>> begin validate tree\n" << std::flush;

	/*
	 * check_size returns T->size on every path, so comparing it against
	 * root_->size could never fail. The real work is the assert inside it,
	 * which is why this only means anything in a debug build.
	 */
	check_size<leaf_type, interior_type>(this->root_);
	std::cout << "Correct size\n" << std::flush;

	// tree property
	if constexpr (is_binary_node<interior_type> ||
		      is_multi_node<interior_type>) {
		if (legal_box(check_box<leaf_type, interior_type>(
			    this->root_, this->tree_box_))) {
			std::cout << "Correct bounding Box\n" << std::flush;
		} else {
			std::cout << "wrong bounding Box\n" << std::flush;
			abort();
		}

		// NOTE: used to check rotate dimension
		// For kdtree binary node, the dummy node may break the rotation
		// manner, since if current dimension is un-splitable, one has
		// to switch to another dimension
		if constexpr (is_rotate_dim_split<
				      typename SplitRule::dim_rule_type> &&
			      is_multi_node<interior_type>) {
			check_tree_same_sequential<leaf_type, interior_type>(
				this->root_, 0);
			std::cout << "Correct rotate dimension\n" << std::flush;
		}
	} else if constexpr (is_dynamic_node<interior_type>) {
		// BUG: below is wrong, the leaf cover circle should centered at
		// one node. Cannot compute it, but to retrieve from the tree
		auto root_cover_circle =
			this->root_->is_leaf
				? get_circle<
					  typename interior_type::CircleType>(
					  static_cast<leaf_type *>(this->root_)
						  ->pts.cut(0,
							    this->root_->size))
				: static_cast<interior_type *>(this->root_)
					  ->get_cover_circle();

		static_assert(std::same_as<decltype(root_cover_circle),
					   typename interior_type::CircleType>);

		if (parlay::all_of(check_cover<leaf_type, interior_type>(
					   this->root_, root_cover_circle),
				   [&](auto const &p) {
					   return within_circle(
						   p, root_cover_circle);
				   })) {
			std::cout << "Correct cover circle\n" << std::flush;
		} else {
			std::cout << "wrong bounding Box\n" << std::flush;
			abort();
		}
	} else {
		assert(false);
	}

	std::cout << "<<< end validate tree\n" << std::flush;
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
size_t base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_tree_height()
{
	if (this->root_ == nullptr) {
		return 0;
	}
	size_t deep = 0;
	return get_max_tree_depth<leaf_type, interior_type>(this->root_, deep);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
size_t base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_max_tree_depth(
	node *T, size_t deep)
{
	if (T == nullptr) {
		return deep;
	}
	if (T->is_leaf) {
		return deep;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		int l = get_max_tree_depth<leaf_type, interior_type>(ti->left,
								     deep + 1);
		int r = get_max_tree_depth<leaf_type, interior_type>(ti->right,
								     deep + 1);
		return std::max(l, r);
	} else if constexpr (is_multi_node<interior_type>) {
		size_t max_depth = 0;
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			max_depth = std::max(
				max_depth,
				get_max_tree_depth<leaf_type, interior_type>(
					ti->tree_nodes[i], deep + 1));
		}
		return max_depth;
	} else {
		return 0;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
double base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_ave_tree_height()
{
	if (this->root_ == nullptr) {
		return 0.0;
	}
	size_t leaf_num = 0, depth_sum = 0;
	count_tree_heights<leaf_type, interior_type>(this->root_, 0, leaf_num,
						     depth_sum);
	return leaf_num == 0 ? 0.0 : double(depth_sum) / double(leaf_num);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
size_t base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::count_tree_nodes_num(
	node *T)
{
	if (T == nullptr) {
		return 0;
	}
	if (T->is_leaf) {
		return 1;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		size_t l, r;
		parlay::par_do(
			[&]() {
				l = count_tree_nodes_num<leaf_type,
							 interior_type>(
					ti->left);
			},
			[&]() {
				r = count_tree_nodes_num<leaf_type,
							 interior_type>(
					ti->right);
			});
		return l + r + 1;
	} else {
		size_t sum = 0;
		for (int i = 0; i < ti->tree_nodes.size(); i++) {
			sum += count_tree_nodes_num<leaf_type, interior_type>(
				ti->tree_nodes[i]);
		}
		return sum + 1;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::count_tree_heights(
	node *T, size_t deep, size_t &leaf_num, size_t &depth_sum)
{
	if (T->is_leaf) {
		leaf_num++;
		depth_sum += deep;
		return;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		count_tree_heights<leaf_type, interior_type>(
			ti->left, deep + 1, leaf_num, depth_sum);
		count_tree_heights<leaf_type, interior_type>(
			ti->right, deep + 1, leaf_num, depth_sum);
	} else if constexpr (is_multi_node<interior_type>) {
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			count_tree_heights<leaf_type, interior_type>(
				ti->tree_nodes[i], deep + 1, leaf_num,
				depth_sum);
		}
	}
}

} // namespace psi

#endif // PSI_BASE_TREE_IMPL_VALIDATION_HPP_
