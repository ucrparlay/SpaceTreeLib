#pragma once
#include <parlay/slice.h>
#include <cstdint>
#include "basic_point.h"

namespace cpdd {
struct alloc_normal_leaf_tag {};
struct alloc_dummy_leaf_tag {};
struct alloc_fat_leaf_tag {};
struct alloc_thin_leaf_tag {};

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
  leaf(slice In, const auto alloc_size, alloc_normal_leaf_tag) :
      node{true, false, static_cast<size_t>(In.size())},
      pts(points::uninitialized(alloc_size)) {
    assert(In.size() <= alloc_size);
    for (int i = 0; i < In.size(); i++) {
      // parlay::assign_dispatch(pts[i], *(In[i].second),
      //                         parlay::uninitialized_copy_tag());
      // pts[i] = In[i].second;
      pts[i] = *(In[i].second);
    }
  }
  // WARN: fat leaf should ensure alloc using parallel copy
  leaf(slice In, alloc_fat_leaf_tag) :
      node{true, false, static_cast<size_t>(In.size())},
      pts(In.begin(), In.end()) {}
  leaf(slice In, alloc_dummy_leaf_tag) :
      node{true, true, static_cast<size_t>(In.size())},
      pts(In.begin(), In.end()) {}
};

// NOTE: point: how data is stored in the tree
// NOTE: slice: input slice
template<typename point, typename slice, typename leaf_alloc_tag>
static leaf<point, slice>* alloc_leaf_node(slice In, const auto alloc_size) {
  using leaf = leaf<point, slice>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new (o) leaf(In, alloc_size, leaf_alloc_tag());
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
