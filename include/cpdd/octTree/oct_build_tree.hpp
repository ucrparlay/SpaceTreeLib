#pragma once

#include <cassert>
#include "../oct_tree.h"

namespace cpdd {

template<typename point>
void
octTree<point>::build( slice A, const dim_type DIM ) {
  points B = points::uninitialized( A.size() );
  this->bbox = this->get_box( A );
  this->root = build_recursive( A, B.cut( 0, A.size() ), 0, DIM, this->bbox );
  assert( this->root != nullptr );
  return;
}

template<typename point>
typename octTree<point>::node*
octTree<point>::serial_build_recursive( slice In, slice Out, dim_type dim,
                                        const dim_type DIM, const box& bx ) {
  node* T;
  return T;
};

template<typename point>
typename octTree<point>::node*
octTree<point>::build_recursive( slice In, slice Out, dim_type dim, const dim_type DIM,
                                 const box& bx ) {
  node* T;
  return T;
};

}  // namespace cpdd
