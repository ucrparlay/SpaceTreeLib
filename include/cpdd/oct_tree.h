#pragma once

#include "base_tree.h"
#include "libmorton/morton.h"

namespace cpdd {

template<typename point>
class octTree : public baseTree<point> {
 public:
  using baseTree = baseTree<point>;
  using node = baseTree::node;
  using leaf = baseTree::leaf;

  using bucket_type = baseTree::bucket_type;
  using balls_type = baseTree::balls_type;
  using dim_type = baseTree::dim_type;
  using z_bit_type = uint_fast8_t;
  using z_value_type = uint_fast64_t;

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

  static inline void interleave_bits(point* p, const dim_type DIM);

  void build(slice In, const dim_type DIM) override;

  void delete_tree() override;

  node* serial_build_recursive(slice In, z_bit_type bit, const dim_type DIM);

  node* build_recursive(slice In, z_bit_type bit, const dim_type DIM);
};

}  // namespace cpdd

#include "octTree/oct_build_tree.hpp"
#include "octTree/oct_inter_node.hpp"
#include "octTree/oct_override.hpp"
