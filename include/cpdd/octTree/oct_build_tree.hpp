#pragma once

#include <cassert>
#include "../oct_tree.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"

namespace cpdd {
template<typename point>
void
octTree<point>::interleave_bits(point* p, const dim_type DIM) {
  z_bit_type loc = 0;
  p->id = 0;
  for (int i = 0; i < KEY_BITS / DIM; i++) {
    for (int d = 0; d < DIM; d++) {
      p->id = p->id |
              (((p->pnt[d] >> i) & static_cast<z_value_type>(1)) << (loc++));
    }
  }
  return;
}

template<typename point>
void
octTree<point>::build(slice A, const dim_type DIM) {
  this->bbox = this->get_box(A);
  parlay::parallel_for(0, A.size(),
                       [&](size_t i) { interleave_bits(&A[i], DIM); });
  parlay::internal::integer_sort_inplace(A,
                                         [&](const point& p) { return p.id; });
  this->root = build_recursive(A, DIM * (KEY_BITS / DIM));
  assert(this->root != nullptr);
  return;
}

template<typename point>
typename octTree<point>::node*
octTree<point>::serial_build_recursive(slice In, z_bit_type bit) {
  node* T;
  return T;
}

template<typename point>
typename octTree<point>::node*
octTree<point>::build_recursive(slice In, z_bit_type bit) {
  node* T;
  return T;
}

}  // namespace cpdd
