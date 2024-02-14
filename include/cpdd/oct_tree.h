#pragma once

#include "comparator.h"
#include "baseTree.h"
#include "basic_point.h"
#include "query_op/nn_search_helpers.h"

namespace cpdd {

#define LOG  std::cout
#define ENDL std::endl << std::flush

template<typename point>
class octTree : public baseTree<point> {
 public:
  using baseTree = baseTree<point>;
  using node = baseTree::node;
  using leaf = baseTree::leaf;

  using bucket_type = baseTree::bucket_type;
  using balls_type = baseTree::balls_type;
  using dim_type = baseTree::dim_type;

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

  virtual void build( slice In, const dim_type DIM ) override;
  virtual node* serial_build_recursive( slice In, slice Out, dim_type dim,
                                        const dim_type DIM, const box& bx ) override;
  virtual node* build_recursive( slice In, slice Out, dim_type dim, const dim_type DIM,
                                 const box& bx ) override;
};

}  // namespace cpdd
