#pragma once

#include <parlay/slice.h>
#include <utility>
#include <algorithm>
#include "../oct_tree.h"
#include "libmorton/morton_BMI.h"
#include "parlay/primitives.h"

namespace cpdd {
template<typename point>
void octTree<point>::build(slice A, const dim_type DIM) {
  // build_z_value(A, DIM);
  // build_z_value_pointer(A, DIM);
  // build_point_z_value(A, DIM);
  build_point(A, DIM);
  return;
}

template<typename point>
void octTree<point>::interleave_bits(point* p, const dim_type DIM) {
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
typename octTree<point>::z_value_type octTree<point>::get_z_value(
    const point& p) {
  assert(point::DIM == 3);
  return libmorton::m3D_e_BMI<z_value_type, coord>(p.pnt[0], p.pnt[1],
                                                   p.pnt[2]);
}

// NOTE: 1: z value only
template<typename point>
void octTree<point>::build_z_value(slice A, const dim_type DIM) {
  parlay::internal::timer t;
  LOG << "here" << ENDL;
  t.start();
  this->bbox = this->get_box(A);
  t.next("get_box");

  auto z_value_arr =
      parlay::map(A, [&](const point& p) { return this->get_z_value(p); });
  t.next("generate_z_value");

  parlay::internal::integer_sort_inplace(
      parlay::make_slice(z_value_arr),
      [&](const z_value_type& val) { return val; });
  t.next("integer_sort_inplace");

  this->root = build_recursive_with_z_value(parlay::make_slice(z_value_arr),
                                            DIM * (KEY_BITS / DIM), DIM);
  t.next("build_recursive_z_value");

  assert(this->root != nullptr);
  return;
}

template<typename point>
node* octTree<point>::build_recursive_with_z_value(z_value_slice In,
                                                   z_bit_type bit,
                                                   const dim_type DIM) {
  size_t n = In.size();
  if (bit == 0 || n < this->LEAVE_WRAP) {
    return alloc_leaf_node<z_value_type, z_value_slice>(
        In, std::max(In.size(), static_cast<size_t>(this->LEAVE_WRAP)));
  }

  z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
  z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
                                  : ~(~(static_cast<z_value_type>(0)) << bit);
  auto less = [&](const z_value_type& x) { return (x & mask) < val; };
  size_t pos = parlay::internal::binary_search(In, less);

  if (pos == 0 || pos == n) {
    return build_recursive_with_z_value(In, bit - 1, DIM);
  }

  node *L, *R;
  parlay::par_do_if(
      n >= this->SERIAL_BUILD_CUTOFF,
      [&] { L = build_recursive_with_z_value(In.cut(0, pos), bit - 1, DIM); },
      [&] { R = build_recursive_with_z_value(In.cut(pos, n), bit - 1, DIM); });
  return alloc_oct_interior_node(L, R, bit);
}

// NOTE: 2: z value and the pointer
template<typename point>
void octTree<point>::build_z_value_pointer(slice A, const dim_type DIM) {
  parlay::internal::timer t;
  LOG << "here" << ENDL;
  t.start();
  this->bbox = this->get_box(A);
  t.next("get_box");

  auto z_value_arr = parlay::map(
      A, [&](point& p) { return std::make_pair(this->get_z_value(p), &p); });
  t.next("generate z_value_pointer");

  parlay::internal::integer_sort_inplace(
      parlay::make_slice(z_value_arr),
      [&](const auto& val) { return val.first; });
  t.next("integer_sort_inplace");

  this->root = build_recursive_with_z_value_pointer(
      parlay::make_slice(z_value_arr), DIM * (KEY_BITS / DIM), DIM);
  t.next("build_recursive_z_value");

  assert(this->root != nullptr);
  return;
}

template<typename point>
node* octTree<point>::build_recursive_with_z_value_pointer(
    z_value_pointer_slice In, z_bit_type bit, const dim_type DIM) {
  size_t n = In.size();
  if (bit == 0 || n < this->LEAVE_WRAP) {
    return alloc_leaf_node<point*, z_value_pointer_slice>(
        In, std::max(In.size(), static_cast<size_t>(this->LEAVE_WRAP)));
  }

  z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
  z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
                                  : ~(~(static_cast<z_value_type>(0)) << bit);
  auto less = [&](const z_value_pointer_pair& x) {
    return (x.first & mask) < val;
  };
  size_t pos = parlay::internal::binary_search(In, less);

  if (pos == 0 || pos == n) {
    return build_recursive_with_z_value_pointer(In, bit - 1, DIM);
  }

  node *L, *R;
  parlay::par_do_if(
      n >= this->SERIAL_BUILD_CUTOFF,
      [&] {
        L = build_recursive_with_z_value_pointer(In.cut(0, pos), bit - 1, DIM);
      },
      [&] {
        R = build_recursive_with_z_value_pointer(In.cut(pos, n), bit - 1, DIM);
      });
  return alloc_oct_interior_node(L, R, bit);
}

// NOTE: 3: sort with z-value
template<typename point>
void octTree<point>::build_point_z_value(slice A, const dim_type DIM) {
  parlay::internal::timer t;
  LOG << "here" << ENDL;
  t.start();
  this->bbox = this->get_box(A);
  t.next("get_box");

  auto z_value_arr = parlay::map(
      A, [&](point& p) { return std::make_pair(this->get_z_value(p), p); });
  t.next("generate z_value_pointer");

  parlay::sort_inplace(parlay::make_slice(z_value_arr));
  t.next("integer_sort_inplace");

  this->root = build_recursive_point_z_value(parlay::make_slice(z_value_arr),
                                             DIM * (KEY_BITS / DIM), DIM);
  t.next("build_recursive");
  assert(this->root != nullptr);
  return;
}

template<typename point>
node* octTree<point>::build_recursive_point_z_value(z_value_point_slice In,
                                                    z_bit_type bit,
                                                    const dim_type DIM) {
  size_t n = In.size();
  if (bit == 0) {
    return alloc_leaf_node<point, z_value_point_slice>(
        In, std::max(In.size(), static_cast<size_t>(this->LEAVE_WRAP)));
  }

  z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
  z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
                                  : ~(~(static_cast<z_value_type>(0)) << bit);
  auto less = [&](const z_value_point_pair& x) {
    return (x.first & mask) < val;
  };
  size_t pos = parlay::internal::binary_search(In, less);

  if (pos == 0 || pos == n) {
    return build_recursive_point_z_value(In, bit - 1, DIM);
  }

  node *L, *R;
  parlay::par_do_if(
      n >= this->SERIAL_BUILD_CUTOFF,
      [&] { L = build_recursive_point_z_value(In.cut(0, pos), bit - 1, DIM); },
      [&] { R = build_recursive_point_z_value(In.cut(pos, n), bit - 1, DIM); });
  return alloc_oct_interior_node(L, R, bit);
}

// NOTE: 4: sort without explicitly computing the z-value
template<typename point>
void octTree<point>::build_point(slice A, const dim_type DIM) {
  parlay::internal::timer t;
  LOG << "here" << ENDL;
  t.start();
  this->bbox = this->get_box(A);
  t.next("get_box");

  t.next("generate z_value_pointer");
  // compares the interleaved bits of points p and q without explicitly
  // interleaving them.  From Timothy Chan.
  auto less = [&](const point& p, const point& q) {
    dim_type j, k;
    coord y, x = 0;
    auto less_msb = [](coord x, coord y) { return x < y && x < (x ^ y); };
    for (j = k = 0; k < DIM; ++k)
      if (less_msb(x, y = p.pnt[k] ^ q.pnt[k])) {
        j = k;
        x = y;
      }
    return p.pnt[j] < q.pnt[j];
  };
  parlay::sort_inplace(A, less);
  t.next("integer_sort_inplace");

  this->root = build_recursive_point(A, DIM * sizeof(coord) * 8 - 1, DIM);
  t.next("build_recursive");
  assert(this->root != nullptr);
  return;
}

template<typename point>
node* octTree<point>::build_recursive_point(slice In, z_bit_type bit,
                                            const dim_type DIM) {
  size_t n = In.size();
  if (bit == 0) {
    return alloc_leaf_node<point, slice>(
        In, std::max(In.size(), static_cast<size_t>(this->LEAVE_WRAP)));
  }

  auto bits = parlay::delayed::map(In, [&](const point& p) {
    return 1 == ((p.pnt[DIM - bit % DIM - 1] >> bit / DIM) & 1);
  });
  size_t pos = std::lower_bound(bits.begin(), bits.end(), 1) - bits.begin();

  if (pos == 0 || pos == n) {
    return build_recursive_point(In, bit - 1, DIM);
  }

  node *L, *R;
  parlay::par_do_if(
      n >= this->SERIAL_BUILD_CUTOFF,
      [&] { L = build_recursive_point(In.cut(0, pos), bit - 1, DIM); },
      [&] { R = build_recursive_point(In.cut(pos, n), bit - 1, DIM); });
  return alloc_oct_interior_node(L, R, bit);
}

}  // namespace cpdd
