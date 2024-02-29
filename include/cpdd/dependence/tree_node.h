#pragma once
#include <parlay/slice.h>
#include <cstdint>
#include "basic_point.h"

namespace cpdd {

struct node {
  node() : is_leaf{false}, is_dummy{false}, size{0} {};

  node(bool _is_leaf, bool _is_dummy, size_t _size) :
      is_leaf{_is_leaf}, is_dummy{_is_dummy}, size{_size} {};

  bool is_leaf;
  bool is_dummy;
  size_t size;
};

template<typename point, typename slice>
struct leaf : node {
  using points = parlay::sequence<point>;
  points pts;
  leaf() : node{true, false, static_cast<size_t>(0)} {};
  // TODO: design better inference in order to distinguish between dummy leaf
  // allocate and size
  leaf(slice In, const size_t size, bool _is_dummy) :
      node{true, false, static_cast<size_t>(In.size())},
      pts(points::uninitialized(size)) {
    assert(In.size() <= size);
    for (int i = 0; i < In.size(); i++) {
      // pts[i] = *(In[i].second);
      pts[i] = In[i].second;
    }
  }
  leaf(slice In, bool _is_dummy) :
      node{true, true, static_cast<size_t>(In.size())},
      pts(points::uninitialized(1)) {
    pts[0] = In[0];
  }
};

// NOTE: point: how data is stored in the tree
// NOTE: slice: input slice
template<typename point, typename slice>
static leaf<point, slice>* alloc_leaf_node(slice In, const size_t alloc_size) {
  using leaf = leaf<point, slice>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new (o) leaf(In, alloc_size, false);
  assert(o->is_dummy == false);
  return o;
}

// template<typename point, typename slice, typename assign_func>
// static leaf<point, slice, assign_func>* alloc_dummy_leaf(
//     slice In,
//     assign_func const& ASSIGN_FUNC = [](point& p, point& q) { p = q; }) {
//   using leaf = leaf<point, slice, assign_func>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf(In, true, ASSIGN_FUNC);
//   assert(o->is_dummy == true);
//   return o;
// }

// template<typename point, typename slice, typename assign_function>
// static leaf<point, slice, assign_function>* alloc_empty_leaf() {
//   using leaf = leaf<point, slice, assign_function>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf();
//   assert(o->size == 0 && o->pts.size() == 0);
//   return o;
// }

// template<typename point>
// static void
// free_leaf(node* T) {
//   parlay::type_allocator<leaf<point>>::retire(static_cast<leaf<point>*>(T));
// }

template<typename point, typename node_type>
static void free_node(node* T) {
  // TODO: add static type check
  parlay::type_allocator<node_type>::retire(static_cast<node_type*>(T));
}

}  // namespace cpdd
