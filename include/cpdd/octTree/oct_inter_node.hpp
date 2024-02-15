#pragma once

#include "../oct_tree.h"

namespace cpdd {
template<typename point>
struct octTree<point>::interior : node {
  node* left;
  node* right;
  interior( node* _left, node* _right ) :
      node{ false, false, _left->size + _right->size }, left( _left ), right( _right ) {}
};

template<typename point>
typename octTree<point>::interior*
octTree<point>::alloc_oct_interior_node( node* L, node* R ) {
  interior* o = parlay::type_allocator<interior>::alloc();
  new ( o ) interior( L, R );
  return o;
}

}  // namespace cpdd
