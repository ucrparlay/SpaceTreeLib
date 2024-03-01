#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"
#include "libmorton/morton.h"

namespace cpdd {

template<typename point>
class octTree : public baseTree<point> {
 public:
  using baseTree = baseTree<point>;

  using bucket_type = baseTree::bucket_type;
  using balls_type = baseTree::balls_type;
  using dim_type = baseTree::dim_type;
  using z_bit_type = uint_fast8_t;
  using z_value_type = uint_fast64_t;
  using z_value_slice = parlay::slice<z_value_type*, z_value_type*>;
  using z_value_pointer_pair = std::pair<z_value_type, point*>;
  using z_value_pointer_slice =
      parlay::slice<z_value_pointer_pair*, z_value_pointer_pair*>;
  using z_value_point_pair = std::pair<z_value_type, point>;
  using z_value_point_slice =
      parlay::slice<z_value_point_pair*, z_value_point_pair*>;

  using coord = typename point::coord;
  using coords = typename point::coords;
  using Num = Num_Comparator<coord>;
  using slice = baseTree::slice;
  using points = baseTree::points;
  using points_iter = baseTree::points_iter;
  using splitter = baseTree::splitter;
  using splitter_s = baseTree::splitter_s;
  using box = baseTree::box;
  using box_s = baseTree::box_s;
  using circle = baseTree::circle;

  static constexpr z_bit_type KEY_BITS = 64;

  struct interior;

  static interior* alloc_oct_interior_node(node* L, node* R, z_bit_type bit);

  static inline z_value_type interleave_bits(point* p, const dim_type DIM);

  static inline z_value_type get_z_value(const point& p);

  void build_z_value(slice In, const dim_type DIM);
  node* build_recursive_with_z_value(z_value_slice In, z_bit_type bit,
                                     const dim_type DIM);

  void build_z_value_pointer(slice In, const dim_type DIM);
  node* build_recursive_with_z_value_pointer(z_value_pointer_slice In,
                                             z_bit_type bit,
                                             const dim_type DIM);

  void build_point_z_value(slice In, const dim_type DIM);
  node* build_recursive_point_z_value(z_value_point_slice In, z_bit_type bit,
                                      const dim_type DIM);

  void build_point(slice In, const dim_type DIM);
  node* build_recursive_point(slice In, z_bit_type bit, const dim_type DIM);

  // NOTE: wrapper
  void build(slice In, const dim_type DIM) override;
  node* serial_build_recursive(slice In, z_bit_type bit, const dim_type DIM);
  node* build_recursive(slice In, z_bit_type bit, const dim_type DIM);

  void delete_tree() override;

  uint64_t binary_search_time = 0;
  uint64_t leaf_alloc_time = 0;
  uint64_t leaf_num = 0;
  uint64_t time_base = 1000000;
};

}  // namespace cpdd

#include "octTree/oct_build_tree.hpp"
#include "octTree/oct_inter_node.hpp"
#include "octTree/oct_override.hpp"
