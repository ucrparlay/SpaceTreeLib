#pragma once
#include <parlay/slice.h>
#include <cstdint>
#include "basic_point.h"
#include "utility.h"

namespace cpdd {

struct node {
    node() : is_leaf{false}, size{0} {};
    node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

    bool is_leaf;
    size_t size;
};

template<typename Point, typename slice, typename point_assign_tag>
struct leaf : node {
    using Points = parlay::sequence<Point>;

    // NOTE: default allocator
    leaf() : node{true, static_cast<size_t>(0)}, is_dummy(false){};

    // NOTE: alloc a normal leaf
    leaf(slice In, const auto alloc_size, alloc_normal_leaf_tag) :
        node{true, static_cast<size_t>(In.size())}, is_dummy(false), pts(Points::uninitialized(alloc_size)) {
        assert(In.size() <= alloc_size);
        if constexpr (is_pair<typename slice::value_type>::value) {  // NOTE: if input is pair, assign using the second
                                                                     // value
            static_assert(std::is_pointer<typename slice::value_type::second_type>::value);
            for (int i = 0; i < In.size(); i++) {
                parlay::assign_dispatch(pts[i], *(In[i].second), point_assign_tag());
            }
        } else {  //  NOTE: input is a Point with coordinates
            for (int i = 0; i < In.size(); i++) {
                parlay::assign_dispatch(pts[i], In[i], point_assign_tag());
            }
        }
    }

    // NOTE: alloc a dummy leaf
    leaf(slice In, alloc_dummy_leaf_tag) :
        node{true, static_cast<size_t>(In.size())}, is_dummy(true), pts(Points::uninitialized(1)) {
        assert(In.size() == 1);
        if constexpr (is_pair<typename slice::value_type>::value) {
            static_assert(std::is_pointer<typename slice::value_type::second_type>::value);
            parlay::assign_dispatch(pts[0], *(In[0].second), point_assign_tag());
        } else {
            parlay::assign_dispatch(pts[0], In[0], point_assign_tag());
        }
    }

    bool is_dummy;
    Points pts;
};

// NOTE: Point: how data is stored in the tree
// NOTE: slice: input slice
template<typename Point, typename slice, typename leaf_alloc_tag, typename point_assign_tag>
static leaf<Point, slice, point_assign_tag>* alloc_leaf_node(slice In, const auto alloc_size) {
    using leaf = leaf<Point, slice, point_assign_tag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, alloc_size, leaf_alloc_tag());
    assert(o->is_dummy == false);
    return o;
}

template<typename Point, typename slice, typename leaf_alloc_tag, typename point_assign_tag>
static leaf<Point, slice, point_assign_tag>* alloc_leaf_node(slice In) {
    using leaf = leaf<Point, slice, point_assign_tag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, leaf_alloc_tag());
    assert(o->is_dummy == false);
    return o;
}

// template<typename Point, typename slice, typename assign_func>
// static leaf<Point, slice, assign_func>* alloc_dummy_leaf(
//     slice In,
//     assign_func const& ASSIGN_FUNC = [](Point& p, Point& q) { p = q; }) {
//   using leaf = leaf<Point, slice, assign_func>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf(In, true, ASSIGN_FUNC);
//   assert(o->is_dummy == true);
//   return o;
// }

// template<typename Point, typename slice, typename assign_function>
// static leaf<Point, slice, assign_function>* alloc_empty_leaf() {
//   using leaf = leaf<Point, slice, assign_function>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf();
//   assert(o->size == 0 && o->pts.size() == 0);
//   return o;
// }

// template<typename Point>
// static void
// free_leaf(node* T) {
//   parlay::type_allocator<leaf<Point>>::retire(static_cast<leaf<Point>*>(T));
// }

template<typename Point, typename node_type>
static void free_node(node* T) {
    // TODO: add static type check
    parlay::type_allocator<node_type>::retire(static_cast<node_type*>(T));
}

}  // namespace cpdd
