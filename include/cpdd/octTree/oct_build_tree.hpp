#pragma once

#include "../oct_tree.h"

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
  parlay::internal::timer t;
  LOG << "here" << ENDL;
  t.start();
  this->bbox = this->get_box(A);
  t.next("get_box");
  parlay::parallel_for(0, A.size(),
                       [&](size_t i) { interleave_bits(&A[i], DIM); });
  t.next("interleave_bits");
  parlay::internal::integer_sort_inplace(A,
                                         [&](const point& p) { return p.id; });
  t.next("integer_sort_inplace");
  this->root = build_recursive(A, DIM * (KEY_BITS / DIM));
  t.next("build_recursive");
  assert(this->root != nullptr);
  return;
}

template<typename point>
typename octTree<point>::node*
octTree<point>::serial_build_recursive(slice In, z_bit_type bit) {
  size_t n = In.size();
  if (bit == 0 || n < this->LEAVE_WRAP) {
    return alloc_leaf_node<point>(In, this->LEAVE_WRAP);
  }

  z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
  z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
                                  : ~(~(static_cast<z_value_type>(0)) << bit);
  auto less = [&](const point& x) { return (x.id & mask) < val; };
  size_t pos = parlay::internal::binary_search(In, less);

  if (pos == 0 || pos == In.size()) {
    return serial_build_recursive(In, bit - 1);
  }

  node *L, *R;
  L = serial_build_recursive(In.cut(0, pos), bit - 1);
  R = serial_build_recursive(In.cut(pos, n), bit - 1);
  return alloc_oct_interior_node(L, R, bit);
}

template<typename point>
typename octTree<point>::node*
octTree<point>::build_recursive(slice In, z_bit_type bit) {
  size_t n = In.size();
  if (bit == 0) {
    return alloc_leaf_node<point>(In, this->LEAVE_WRAP);
  }
  if (n < this->SERIAL_BUILD_CUTOFF) {
    return serial_build_recursive(In, bit);
  }

  z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
  z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
                                  : ~(~(static_cast<z_value_type>(0)) << bit);
  auto less = [&](const point& x) { return (x.id & mask) < val; };
  size_t pos = parlay::internal::binary_search(In, less);

  if (pos == 0 || pos == n) {
    return build_recursive(In, bit - 1);
  }

  node *L, *R;
  parlay::par_do([&] { L = build_recursive(In.cut(0, pos), bit - 1); },
                 [&] { R = build_recursive(In.cut(pos, n), bit - 1); });
  return alloc_oct_interior_node(L, R, bit);
}

}  // namespace cpdd
