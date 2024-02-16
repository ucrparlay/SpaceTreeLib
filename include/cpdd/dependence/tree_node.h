#pragma once
#include <cstdint>
#include "basic_point.h"

namespace cpdd {

template<typename point>
struct node {
  using coord = point::coord;
  using coords = point::coords;

  node() : is_leaf{false}, is_dummy{false}, size{0} {};

  node(bool _is_leaf, bool _is_dummy, size_t _size) :
      is_leaf{_is_leaf}, is_dummy{_is_dummy}, size{_size} {};

  bool is_leaf;
  bool is_dummy;
  size_t size;
};

template<typename point>
struct leaf : node<point> {
  using node = node<point>;
  using points = parlay::sequence<point>;
  using slice = parlay::slice<point*, point*>;
  points pts;
  leaf() : node{true, false, static_cast<size_t>(0)} {};
  leaf(slice In, const size_t size) :
      node{true, false, static_cast<size_t>(In.size())} {
    assert(In.size() <= size);
    pts = points::uninitialized(size);
    for (int i = 0; i < In.size(); i++) {
      pts[i] = In[i];
    }
  }
  leaf(slice In, bool _is_dummy) :
      node{true, true, static_cast<size_t>(In.size())} {
    pts = points::uninitialized(1);
    pts[0] = In[0];
  }
};

template<typename point, typename slice = parlay::slice<point*, point*>>
static leaf<point>*
alloc_leaf_node(slice In, const uint_fast8_t LEAVE_WRAP) {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new (o) leaf(In, In.size() <= LEAVE_WRAP ? LEAVE_WRAP : In.size());
  assert(o->is_dummy == false);
  return o;
}

template<typename point, typename slice = parlay::slice<point*, point*>>
static leaf<point>*
alloc_dummy_leaf(slice In) {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new (o) leaf(In, true);
  assert(o->is_dummy == true);
  return o;
}

template<typename point, typename slice = parlay::slice<point*, point*>>
static leaf<point>*
alloc_empty_leaf() {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new (o) leaf();
  assert(o->size == 0 && o->pts.size() == 0);
  return o;
}

template<typename point>
static void
free_leaf(node<point>* T) {
  parlay::type_allocator<leaf<point>>::retire(static_cast<leaf<point>*>(T));
}

template<typename point, typename node_type>
static void
free_node(node<point>* T) {
  // TODO: add static type check
  parlay::type_allocator<node_type>::retire(static_cast<node_type*>(T));
}

}  // namespace cpdd
