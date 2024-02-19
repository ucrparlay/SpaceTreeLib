#pragma once

#include <parlay/slice.h>
#include "../oct_tree.h"
#include "libmorton/morton_BMI.h"
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
typename octTree<point>::z_value_type
octTree<point>::get_z_value(const point& p) {
  assert(point::DIM == 3);
  return libmorton::m3D_e_BMI<z_value_type, coord>(p.pnt[0], p.pnt[1],
                                                   p.pnt[2]);
}

template<typename point>
void
octTree<point>::build(slice A, const dim_type DIM) {
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
node*
octTree<point>::build_recursive_with_z_value(z_slice In, z_bit_type bit,
                                             const dim_type DIM) {
  size_t n = In.size();
  if (bit == 0 || n < this->LEAVE_WRAP) {
    return alloc_leaf_node<z_value_type, z_slice>(In, this->LEAVE_WRAP);
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

// template<typename point>
// void
// octTree<point>::build(slice A, const dim_type DIM) {
//   parlay::internal::timer t;
//   LOG << "here" << ENDL;
//   t.start();
//   this->bbox = this->get_box(A);
//   t.next("get_box");
//
//   // compares the interleaved bits of points p and q without explicitly
//   // interleaving them.  From Timothy Chan.
//   auto less = [&](const point& p, const point& q) {
//     dim_type j, k;
//     coord y, x = 0;
//     auto less_msb = [](coord x, coord y) { return x < y && x < (x ^ y); };
//     for (j = k = 0; k < DIM; ++k)
//       if (less_msb(x, y = p.pnt[k] ^ q.pnt[k])) {
//         j = k;
//         x = y;
//       }
//     return p.pnt[j] < q.pnt[j];
//   };
//   parlay::sort_inplace(A, less);
//   t.next("integer_sort_inplace");
//
//   this->root = build_recursive(A, DIM * sizeof(coord) * 8 - 1, DIM);
//   t.next("build_recursive");
//   assert(this->root != nullptr);
//   return;
// }

// template<typename point>
// node*
// octTree<point>::serial_build_recursive(slice In, z_bit_type bit,
//                                        const dim_type DIM) {
//   size_t n = In.size();
//   if (bit == 0 || n < this->LEAVE_WRAP) {
//     return alloc_leaf_node<point, z_slice>(In, this->LEAVE_WRAP);
//   }
//
//   auto bits = parlay::delayed::map(In, [&](const point& p) {
//     return 1 == ((p.pnt[DIM - bit % DIM - 1] >> bit / DIM) & 1);
//   });
//   size_t pos = std::lower_bound(bits.begin(), bits.end(), 1) - bits.begin();
//   // z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
//   // z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
//   //                                 : ~(~(static_cast<z_value_type>(0)) <<
//   //                                 bit);
//   // auto less = [&](const point& x) { return (x.id & mask) < val; };
//   // size_t pos = parlay::internal::binary_search(In, less);
//
//   if (pos == 0 || pos == In.size()) {
//     return serial_build_recursive(In, bit - 1, DIM);
//   }
//
//   node *L, *R;
//   L = serial_build_recursive(In.cut(0, pos), bit - 1, DIM);
//   R = serial_build_recursive(In.cut(pos, n), bit - 1, DIM);
//   return alloc_oct_interior_node(L, R, bit);
// }

// template<typename point>
// node*
// octTree<point>::build_recursive(slice In, z_bit_type bit, const dim_type DIM)
// {
//   size_t n = In.size();
//   if (bit == 0) {
//     return alloc_leaf_node<point, z_slice>(In, this->LEAVE_WRAP);
//   }
//   if (n < this->SERIAL_BUILD_CUTOFF) {
//     return serial_build_recursive(In, bit, DIM);
//   }
//
//   auto bits = parlay::delayed::map(In, [&](const point& p) {
//     return 1 == ((p.pnt[DIM - bit % DIM - 1] >> bit / DIM) & 1);
//   });
//   size_t pos = std::lower_bound(bits.begin(), bits.end(), 1) - bits.begin();
//   // z_value_type val = (static_cast<z_value_type>(1)) << (bit - 1);
//   // z_value_type mask = (bit == 64) ? ~(static_cast<z_value_type>(0))
//   //                                 : ~(~(static_cast<z_value_type>(0)) <<
//   //                                 bit);
//   // auto less = [&](const point& x) { return (x.id & mask) < val; };
//   // size_t pos = parlay::internal::binary_search(In, less);
//
//   if (pos == 0 || pos == n) {
//     return build_recursive(In, bit - 1, DIM);
//   }
//
//   node *L, *R;
//   parlay::par_do([&] { L = build_recursive(In.cut(0, pos), bit - 1, DIM); },
//                  [&] { R = build_recursive(In.cut(pos, n), bit - 1, DIM); });
//   return alloc_oct_interior_node(L, R, bit);
// }

}  // namespace cpdd
