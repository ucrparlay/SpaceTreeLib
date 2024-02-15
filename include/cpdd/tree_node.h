#pragma once
#include "basic_point.h"

namespace cpdd {

template<typename point>
struct node {
  using coord = point::coord;
  using coords = point::coords;

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
  leaf() : node{ true, false, static_cast<size_t>( 0 ) } {};
  leaf( slice In ) : node{ true, false, static_cast<size_t>( In.size() ) } {
    pts = points::uninitialized( 32 );
    for ( int i = 0; i < In.size(); i++ ) {
      pts[i] = In[i];
    }
  }
  leaf( slice In, bool _is_dummy ) :
      node{ true, true, static_cast<size_t>( In.size() ) } {
    pts = points::uninitialized( 1 );
    pts[0] = In[0];
  }
};

template<typename point, typename slice = parlay::slice<point*, point*>>
leaf<point>*
alloc_leaf_node( slice In ) {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new ( o ) leaf( In );
  assert( o->is_dummy == false );
  return o;
}

template<typename point, typename slice = parlay::slice<point*, point*>>
leaf<point>*
alloc_dummy_leaf( slice In ) {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new ( o ) leaf( In, true );
  assert( o->is_dummy == true );
  return o;
}

template<typename point, typename slice = parlay::slice<point*, point*>>
leaf<point>*
alloc_empty_leaf() {
  using leaf = leaf<point>;
  leaf* o = parlay::type_allocator<leaf>::alloc();
  new ( o ) leaf();
  assert( o->size == 0 && o->pts.size() == 0 );
  return o;
}

// template<typename point>
// typename baseTree<point>::interior*
// baseTree<point>::alloc_interior_node( node* L, node* R, const splitter& split ) {
//   interior* o = parlay::type_allocator<interior>::alloc();
//   new ( o ) interior( L, R, split );
//   return o;
// }
//
// template<typename point>
// typename baseTree<point>::simple_node*
// baseTree<point>::alloc_simple_node( simple_node* L, simple_node* R ) {
//   simple_node* o = parlay::type_allocator<simple_node>::alloc();
//   new ( o ) simple_node( L, R );
//   return o;
// }
//
// template<typename point>
// typename baseTree<point>::simple_node*
// baseTree<point>::alloc_simple_node( size_t sz ) {
//   simple_node* o = parlay::type_allocator<simple_node>::alloc();
//   new ( o ) simple_node( sz );
//   return o;
// }

template<typename point>
void
free_leaf( node<point>* T ) {
  parlay::type_allocator<leaf<point>>::retire( static_cast<leaf<point>*>( T ) );
}

template<typename point, typename node_type>
void
free_node( node<point>* T ) {
  // TODO: add static type check
  parlay::type_allocator<node_type>::retire( static_cast<node_type*>( T ) );
}

}  // namespace cpdd
