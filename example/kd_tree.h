#pragma once

#include <iostream>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"

// Example usage of KdTree from SpaceTreeLib
// This demonstrates basic operations: Build, KNN, RangeQuery, BatchInsert, and
// BatchDelete

namespace kd_tree_example {

using Coord = long;

// Define augmentation structure for points (stores an ID)
// WARN: All functions must be defined
struct AugId {
  using IdType = int;
  IdType id;

  bool operator<(AugId const& rhs) const { return id < rhs.id; }
  bool operator==(AugId const& rhs) const { return id == rhs.id; }
  friend std::ostream& operator<<(std::ostream& os, AugId const& rhs) {
    os << rhs.id;
    return os;
  }
};

// Define point type: 2D points with augmented ID
using Point = psi::AugPoint<Coord, 2, AugId>;
using Points = parlay::sequence<Point>;
using BT = psi::BaseTree<Point>;

// Leaf augmentation: stores bounding box
// WARN: All functions must be defined
template <class BaseTree>
struct LeafAugBox {
  using Box = typename BaseTree::Box;
  using Slice = typename BaseTree::Slice;

  LeafAugBox() : box(BaseTree::GetEmptyBox()) {}
  LeafAugBox(Box const& _box) : box(_box) {}
  LeafAugBox(Slice In) : box(BaseTree::GetBox(In)) {}

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  void UpdateAug(Slice In) { this->box = BaseTree::GetBox(In); }

  void Reset() { this->box = BaseTree::GetEmptyBox(); }

  Box box;
};

// Interior node augmentation: stores bounding box and parallel build flag
// WARN: All functions must be defined
template <class BaseTree>
struct InteriorAugBox {
  using Box = typename BaseTree::Box;

  InteriorAugBox() : box(BaseTree::GetEmptyBox()) {
    force_par_indicator.reset();
  }
  InteriorAugBox(Box const& _box) : box(_box) { force_par_indicator.reset(); }

  // Binary create
  template <typename Leaf, typename Interior>
  static Box Create(psi::Node* l, psi::Node* r) {
    return BaseTree::GetBox(BaseTree::template RetriveBox<Leaf, Interior>(l),
                            BaseTree::template RetriveBox<Leaf, Interior>(r));
  }

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  // Binary update
  template <typename Leaf, typename Interior>
  void Update(psi::Node* l, psi::Node* r) {
    this->box = Create<Leaf, Interior>(l, r);
  }

  void SetParallelFlag(bool const flag) {
    this->force_par_indicator.emplace(flag);
  }

  void ResetParallelFlag() { this->force_par_indicator.reset(); }

  bool GetParallelFlagIniStatus() {
    return this->force_par_indicator.has_value();
  }

  bool ForceParallel(size_t sz) const {
    return this->force_par_indicator.has_value()
               ? this->force_par_indicator.value()
               : sz > BaseTree::kSerialBuildCutoff;
  }

  void Reset() {
    this->force_par_indicator.reset();
    this->box = BaseTree::GetEmptyBox();
  }

  Box box;
  std::optional<bool> force_par_indicator;
};

// Define split rule: max stretch dimension + object median
using SplitRule = psi::OrthogonalSplitRule<psi::MaxStretchDim<Point>,
                                           psi::ObjectMedian<Point>>;

// Alternative split rule: rotate dimension + spatial median
using AnotherSplitRule =
    psi::OrthogonalSplitRule<psi::RotateDim<Point>, psi::SpatialMedian<Point>>;

// Define KdTree type
using Tree = psi::KdTree<Point, SplitRule, LeafAugBox<BT>, InteriorAugBox<BT>>;

void run_example() {
  // 1. Create sample 2D points
  size_t n = 1000;
  Points points(n);

  // Generate random points in a 1000x1000 grid
  parlay::parallel_for(0, n, [&](size_t i) {
    points[i][0] = (i * 7) % 1000;   // x coordinate
    points[i][1] = (i * 13) % 1000;  // y coordinate
    points[i].aug.id = i;            // unique ID
  });

  std::cout << "Created " << n << " random 2D points" << std::endl;

  // 2. Build the KdTree
  Tree tree;
  auto points_copy =
      points;  // WARN: Important! The input array will be changed during build.
  tree.Build(parlay::make_slice(points_copy));
  std::cout << "Built KdTree with " << n << " points" << std::endl;

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
            << " nearest neighbors to point (500, 500):" << std::endl;
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

}  // namespace kd_tree_example
