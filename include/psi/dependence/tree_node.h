#ifndef PSI_DEPENDENCE_TREE_NODE_H_
#define PSI_DEPENDENCE_TREE_NODE_H_

#include <parlay/slice.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <variant>

#include "psi/dependence/basic_point.h"
#include "psi/dependence/comparator.h"
#include "psi/dependence/concepts.h"
#include "parlay/utilities.h"

namespace psi
{

struct alloc_normal_leaf_tag {
};
struct alloc_dummy_leaf_tag {
};
struct alloc_empty_leaf_tag {
};

struct node {
	node() : is_leaf{false}, size{0} {};
	node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

	// Adding a virtual destructor makes node polymorphic
	virtual ~node() = default;

	bool is_leaf;
	size_t size;
};

template <typename Point, typename Range, uint_fast8_t DefaultWrap,
	  typename AugType = std::monostate,
	  typename PointAssignTag = parlay::move_assign_tag>
struct leaf_node : node {
	using points_type = parlay::sequence<Point>;
	using at_type = AugType;

	// NOTE: default allocator
	leaf_node()
	    : node{true, static_cast<size_t>(0)}, is_dummy(false),
	      aug(AugType())
	{
	}

	// NOTE: alloc a leaf with default size
	leaf_node(Range in, alloc_normal_leaf_tag)
	    : node{true, static_cast<size_t>(in.size())}, is_dummy(false),
	      pts(points_type::uninitialized(DefaultWrap)), aug(at_type(in))
	{
		assert(in.size() <= DefaultWrap);
		std::ranges::for_each(in, [&, i = 0](auto &&x) mutable {
			parlay::assign_dispatch(pts[i++], x, PointAssignTag());
		});
	}

	// NOTE: alloc a leaf with specific size
	leaf_node(Range in, size_t const alloc_size, alloc_normal_leaf_tag)
	    : node{true, static_cast<size_t>(in.size())}, is_dummy(false),
	      pts(points_type::uninitialized(alloc_size)), aug(at_type(in))
	{
		assert(in.size() <= alloc_size);
		std::ranges::for_each(in, [&, i = 0](auto &&x) mutable {
			parlay::assign_dispatch(pts[i++], x, PointAssignTag());
		});
	}

	// NOTE: alloc a dummy leaf
	leaf_node(Range in, alloc_dummy_leaf_tag)
	    : node{true, static_cast<size_t>(in.size())}, is_dummy(true),
	      pts(points_type::uninitialized(1)), aug(at_type(in.cut(0, 1)))
	{
		parlay::assign_dispatch(pts[0], in[0], PointAssignTag());
	}

	points_type const get_points() const
	{
		return pts.substr(0, this->size);
	}

	decltype(auto) get_box()
		requires(node_has_non_trivial_aug<at_type> && has_box<at_type>)
	{
		return aug.get_box();
	}

	decltype(auto) get_box() const
		requires(node_has_non_trivial_aug<at_type> && has_box<at_type>)
	{
		return aug.get_box();
	}

	// TODO: should put this one in the aug as well and don't use the
	// monostate
	auto update_aug(Range in)
	{
		return aug.update_aug(in);
	}

	auto reset_aug()
	{
		return aug.reset();
	}

	bool is_dummy;
	points_type pts;
	AugType aug;
};

// NOTE:: Alloc a leaf with input IN and given size
// TODO: the input aug is no longer valid
template <typename Range, typename leaf_type>
static leaf_type *alloc_fix_size_leaf_node(Range in, size_t const alloc_size)
{
	leaf_type *o = parlay::type_allocator<leaf_type>::alloc();
	new (o) leaf_type(in, alloc_size, alloc_normal_leaf_tag());
	assert(o->is_dummy == false);
	return o;
}

// NOTE: Alloc a leaf
template <typename Range, typename leaf_type>
static leaf_type *alloc_normal_leaf_node(Range in)
{
	leaf_type *o = parlay::type_allocator<leaf_type>::alloc();
	new (o) leaf_type(in, alloc_normal_leaf_tag());
	assert(o->is_dummy == false);
	return o;
}

// NOTE: Alloc a empty leaf_type
template <typename Range, typename leaf_type>
static leaf_type *alloc_empty_leaf_node()
{
	leaf_type *o = parlay::type_allocator<leaf_type>::alloc();
	new (o) leaf_type();
	return o;
}

// NOTE: Alloc a dummy leaf
template <typename Range, typename leaf_type>
static leaf_type *alloc_dummy_leaf_node(Range in)
{
	leaf_type *o = parlay::type_allocator<leaf_type>::alloc();
	new (o) leaf_type(in, alloc_dummy_leaf_tag());
	assert(o->is_dummy == true);
	assert(o->pts.size() == 1);
	return o;
}

template <typename Point, typename SplitType, typename AugType>
struct binary_node : node {
	using pt_type = Point;
	using st_type = SplitType;
	using at_type = AugType;

	binary_node(node *_left, node *_right, st_type const &_split,
		    at_type const &_aug)
	    : node{false, _left->size + _right->size}, left(_left),
	      right(_right), split(_split), aug(_aug)
	{
	}

	// Adding a virtual destructor makes node polymorphic
	virtual ~binary_node() override = default;

	// NOTE: test whether we can fetch @depth levels from @T
	template <typename interior_type>
	static inline bool test_depth(node *T, int cur_depth, int depth)
	{
		if (cur_depth == depth) {
			return true;
		}
		if (T->is_leaf) {
			return false;
		}
		auto ti = static_cast<interior_type *>(T);
		return test_depth<interior_type>(ti->left, cur_depth + 1,
						 depth) &&
		       test_depth<interior_type>(ti->right, cur_depth + 1,
						 depth);
	}

	st_type const &get_split() const
	{
		return split;
	}

	node *left;
	node *right;
	st_type split;
	at_type aug;
};

// NOTE: multi-way node whose children # is fixed
template <typename Point, uint_fast8_t md, typename SplitType, typename AugType>
struct multi_node : node {
	using bucket_type = uint_fast8_t;
	using coord_type = typename Point::coord_type;
	using num_type = num_comparator<coord_type>;
	using st_type = SplitType;
	using at_type = AugType;

	static consteval auto get_regions()
	{
		return 1 << md;
	}

	static consteval auto get_levels()
	{
		return md;
	}

	static consteval auto get_split_nums()
	{
		return md;
	}

	// NOTE: whether we use same splitter for same level
	// e.g., the horizontal in lavel is not the same
	static consteval auto equal_split()
		requires std::same_as<
			st_type, std::array<typename st_type::value_type, md>>
	{
		return true; // every dimension one splitter
	}

	static consteval auto equal_split()
		requires std::same_as<
			st_type,
			std::array<typename st_type::value_type, 1 << md>>
	{
		return false; // the spliiter number mathes inter nodes num
	}

	using node_arr_type = std::array<node *, get_regions()>;

	multi_node(node_arr_type const &_tree_nodes, st_type const &_split,
		   at_type const &_aug)
	    : node{false,
		   std::accumulate(_tree_nodes.begin(), _tree_nodes.end(),
				   static_cast<size_t>(0),
				   [](size_t acc, node *n) -> size_t {
					   return acc + n->size;
				   })},
	      tree_nodes(_tree_nodes), split(_split), aug(_aug)
	{
	}

	// Adding a virtual destructor makes node polymorphic
	// TODO: check whether it is possible to remove it then knn-mix
	virtual ~multi_node() override = default;

	inline size_t merge_size(bucket_type const idx)
	{
		if (idx == 1) {
			return this->size;
		} else if (idx >= get_regions()) {
			return tree_nodes[idx - get_regions()]->size;
		} else {
			return merge_size(2 * idx) + merge_size(2 * idx + 1);
		}
	}

	st_type const &get_split() const
	{
		return split;
	}

	node_arr_type tree_nodes;
	st_type split;
	at_type aug;
};

// NOTE: dynamic multi-way node whose children # is dynamic
template <typename Point, typename SplitType, typename AugType>
struct dynamic_node : node {
	using bucket_type = uint_fast8_t;
	using coord_type = typename Point::coord_type;
	using num_type = num_comparator<coord_type>;
	using st_type = SplitType;
	using at_type = AugType;
	using node_arr_type = parlay::sequence<node *>;

	dynamic_node()
	    : node{false, static_cast<size_t>(0)}, tree_nodes(node_arr_type()),
	      split(st_type()), aug(at_type())
	{
	}

	dynamic_node(node_arr_type const &_tree_nodes, st_type const &_split,
		     at_type const &_aug)
	    : node{false,
		   std::accumulate(_tree_nodes.begin(), _tree_nodes.end(),
				   static_cast<size_t>(0),
				   [](size_t acc, node *n) -> size_t {
					   return acc + n->size;
				   })},
	      tree_nodes(_tree_nodes), split(_split), aug(_aug)
	{
	}

	// Adding a virtual destructor makes node polymorphic
	// TODO: check whether it is possible to remove it then knn-mix
	virtual ~dynamic_node() override = default;

	st_type const &get_split() const
	{
		return split;
	}

	node_arr_type tree_nodes;
	st_type split;
	at_type aug;
};

template <typename interior_type>
static interior_type *
alloc_interior_node(node *L, node *R,
		    typename interior_type::st_type const &split)
{
	interior_type *o = parlay::type_allocator<interior_type>::alloc();
	new (o) interior_type(L, R, split);
	return o;
}

template <typename interior_type>
static interior_type *
alloc_interior_node(node *L, node *R,
		    typename interior_type::st_type const &split,
		    typename interior_type::at_type const &aug)
{
	interior_type *o = parlay::type_allocator<interior_type>::alloc();
	new (o) interior_type(L, R, split, aug);
	return o;
}

template <typename interior_type>
static interior_type *
alloc_interior_node(typename interior_type::node_arr_type const &tree_nodes,
		    typename interior_type::st_type const &split)
{
	interior_type *o = parlay::type_allocator<interior_type>::alloc();
	new (o) interior_type(tree_nodes, split);
	return o;
}

// TODO: maybe replace the node_arr_type, st_type, and at_type with template
// parameters
template <typename interior_type>
static interior_type *
alloc_interior_node(typename interior_type::node_arr_type const &tree_nodes,
		    typename interior_type::st_type const &split,
		    typename interior_type::at_type const &aug)
{
	interior_type *o = parlay::type_allocator<interior_type>::alloc();
	new (o) interior_type(tree_nodes, split, aug);
	return o;
}

template <typename interior_type>
static interior_type *alloc_interior_node()
{
	interior_type *o = parlay::type_allocator<interior_type>::alloc();
	new (o) interior_type();
	return o;
}

template <typename NodeType>
static void free_node(node *T)
{
	parlay::type_allocator<NodeType>::retire(static_cast<NodeType *>(T));
}

// Define some alias for concept check
template <typename T>
concept is_binary_node = check_binary_node<T, psi::binary_node>;

template <typename T>
concept is_multi_node = check_multi_node<T, psi::multi_node>;

template <typename T>
concept is_dynamic_node = check_dynamic_node<T, psi::dynamic_node>;

// TODO: may remove below from inner tree
template <typename T>
concept is_pointer_to_node = check_pointer_to_node<T, node>;

template <typename T, typename Point>
concept is_node_box = requires {
	requires is_pair<T> && is_pointer_to_node<typename T::first_type> &&
			 is_box<typename T::second_type, Point>;
};

} // namespace psi

#endif // PSI_DEPENDENCE_TREE_NODE_H_
