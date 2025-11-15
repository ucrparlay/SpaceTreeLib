#pragma once

#include <iostream>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/p_tree.h"

// Example usage of PTree from SpaceTreeLib
// PTree uses space-filling curves (Morton/Hilbert) to order points

namespace p_tree_example {

using Coord = long;

// Define augmentation structure for points with space-filling curve code
// WARN: All functions must be defined
// PTree requires AugIdCode which includes both an id and a curve code
struct AugIdCode {
  using IdType = int_fast32_t;
  using CurveCode = uint64_t;

  AugIdCode() : code(0), id(0) {}

  void SetMember(CurveCode const& val) { code = val; }

  bool operator<(AugIdCode const& rhs) const {
    return code == rhs.code ? id < rhs.id : code < rhs.code;
  }

  bool operator==(AugIdCode const& rhs) const {
    // WARN: code is not important, we only need to ensure the id
    return id == rhs.id;
  }

  friend std::ostream& operator<<(std::ostream& os, AugIdCode const& rhs) {
    os << rhs.code << " " << rhs.id;
    return os;
  }

  CurveCode code;
  IdType id;
};

// Define point type: 2D points with augmented ID and curve code
using Point = psi::AugPoint<Coord, 2, AugIdCode>;
using Points = parlay::sequence<Point>;
using BT = psi::BaseTree<Point>;

// Define split rule using space-filling curve (Morton curve in this example)
using SplitRule = psi::SpacialFillingCurve<psi::MortonCurve<Point>>;

// Alternative: HilbertCurve<Point>
using AnotherSplitRule = psi::SpacialFillingCurve<psi::HilbertCurve<Point>>;

// Define PTree type
// PTree doesn't use LeafAug/InteriorAug like KdTree/OrthTree
// It uses CPAM (Compressed Purely Functional Augmented Maps) internally to
// maintain the bounding box as an augmented value.
using Tree = psi::PTree<Point, AnotherSplitRule>;

void run_example() {
  // 1. Create sample 2D points
  size_t n = 1000;
  Points points(n);

  // Generate random points in a 1000x1000 grid
  parlay::parallel_for(0, n, [&](size_t i) {
    points[i][0] = (i * 7) % 1000;   // x coordinate
    points[i][1] = (i * 13) % 1000;  // y coordinate
    points[i].aug.id = i;            // unique ID
    // Note: curve code will be computed during build
  });

  std::cout << "Created " << n << " random 2D points" << std::endl;

  // 2. Build the PTree
  Tree tree;
  auto points_copy =
      points;  // WARN: Important! The input array will be changed during build.
  tree.Build(parlay::make_slice(points_copy));
  std::cout << "Built PTree with " << n << " points using Morton curve"
            << std::endl;

  // 3. K-Nearest Neighbors query
  int K = 10;
  Point query_point;
  query_point[0] = 500;
  query_point[1] = 500;
  query_point.aug.id = -1;

  using DisType = typename Point::DisType;
  using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;

  parlay::sequence<nn_pair> knn_result(K, nn_pair(std::ref(points[0]), 0));
  psi::kBoundedQueue<Point, nn_pair> bq(parlay::make_slice(knn_result));

  auto* root = tree.GetRoot();
  tree.KNN(root, query_point, bq);

  std::cout << "Found " << K
            << " nearest neighbors to point (500, 500) (unsorted):"
            << std::endl;
  for (int i = 0; i < std::min(3, K); i++) {
    auto& [pt, dist] = knn_result[i];
    std::cout << "  Point " << pt.get().aug.id << " at (" << pt.get()[0] << ", "
              << pt.get()[1] << ") with distance " << dist << std::endl;
  }

  // 4. Range query (Range count is available as well)
  typename Tree::Box query_box;
  query_box.first[0] = 400;
  query_box.first[1] = 400;
  query_box.second[0] = 600;
  query_box.second[1] = 600;

  Points range_result(n);  // Allocate max possible size
  auto [count, logger] =
      tree.RangeQuery(query_box, parlay::make_slice(range_result));

  std::cout << "Range query [400,600]x[400,600] found " << count << " points"
            << std::endl;

  // 5. Batch insert new points
  size_t insert_count = 100;
  Points new_points(insert_count);
  parlay::parallel_for(0, insert_count, [&](size_t i) {
    new_points[i][0] = (n + i) * 11 % 1000;
    new_points[i][1] = (n + i) * 17 % 1000;
    new_points[i].aug.id = n + i;
  });

  Points insert_copy(new_points);  // same as build
  tree.BatchInsert(parlay::make_slice(insert_copy));
  std::cout << "Inserted " << insert_count << " new points" << std::endl;

  // 6. Batch delete some points
  size_t delete_count = 50;
  Points points_to_delete = points.subseq(0, delete_count);

  Points delete_copy(points_to_delete);  // same as build
  tree.BatchDelete(parlay::make_slice(delete_copy));
  std::cout << "Deleted " << delete_count << " points" << std::endl;

  // 7. Clean up
  tree.DeleteTree();
  std::cout << "Example completed successfully!" << std::endl;
}

}  // namespace p_tree_example
