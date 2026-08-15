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

namespace kd_tree_example
{

// Type for each coordinate
using Coord = long;

// Define augmentation structure for points (stores an ID)
// WARN: All functions must be defined
struct AugId {
	using IdType = int;
	IdType id;

	bool operator<(AugId const &rhs) const
	{
		return id < rhs.id;
	}
	bool operator==(AugId const &rhs) const
	{
		return id == rhs.id;
	}
	friend std::ostream &operator<<(std::ostream &os, AugId const &rhs)
	{
		os << rhs.id;
		return os;
	}
};

// Define point type: 2D points with augmented ID
using Point = psi::AugPoint<Coord, 2, AugId>;
using Points = parlay::sequence<Point>;
using BT = psi::BaseTree<Point>;

// Leaf augmentation: stores bounding box

// Interior node augmentation: stores bounding box and parallel build flag

// Define split rule: max stretch dimension + object median
using SplitRule = psi::OrthogonalSplitRule<psi::MaxStretchDim<Point>,
					   psi::ObjectMedian<Point>>;

// Alternative split rule: rotate dimension + spatial median
using AnotherSplitRule = psi::OrthogonalSplitRule<psi::RotateDim<Point>,
						  psi::SpatialMedian<Point>>;

// Define KdTree type
using Tree = psi::KdTree<Point, SplitRule>;

void run_example()
{
	// 1. Create sample 2D points
	size_t n = 1000;
	Points points(n);

	// Generate random points in a 1000x1000 grid
	parlay::parallel_for(0, n, [&](size_t i) {
		points[i][0] = (i * 7) % 1000;	// x coordinate
		points[i][1] = (i * 13) % 1000; // y coordinate
		points[i].aug.id = i;		// unique ID
	});

	std::cout << "Created " << n << " random 2D points" << std::endl;

	// 2. Build the KdTree
	Tree tree;
	tree.Build(parlay::make_slice(points));
	std::cout << "Built KdTree with " << n << " points" << std::endl;

	// 3. K-Nearest Neighbors query
	int K = 10;
	Point query_point;
	query_point[0] = 500;
	query_point[1] = 500;
	query_point.aug.id = -1;

	using DisType = typename Point::DisType;
	using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;

	parlay::sequence<nn_pair> knn_result(K,
					     nn_pair(std::ref(points[0]), 0));
	psi::kBoundedQueue<Point, nn_pair> bq(parlay::make_slice(knn_result));

	tree.KNN(query_point, bq);

	std::cout << "Found " << K
		  << " nearest neighbors to point (500, 500) (unsorted):"
		  << std::endl;
	for (int i = 0; i < std::min(3, K); i++) {
		auto &[pt, dist] = knn_result[i];
		std::cout << "  Point " << pt.get().aug.id << " at ("
			  << pt.get()[0] << ", " << pt.get()[1]
			  << ") with distance " << dist << std::endl;
	}

	// 4. Range query (Range count is available as well)
	typename Tree::Box query_box;
	query_box.first[0] = 400;
	query_box.first[1] = 400;
	query_box.second[0] = 600;
	query_box.second[1] = 600;

	Points range_result(n); // Allocate max possible size
	auto [count, logger] =
		tree.RangeQuery(query_box, parlay::make_slice(range_result));

	std::cout << "Range query [400,600]x[400,600] found " << count
		  << " points" << std::endl;

	// 5. Batch insert new points
	size_t insert_count = 100;
	Points new_points(insert_count);
	parlay::parallel_for(0, insert_count, [&](size_t i) {
		new_points[i][0] = (n + i) * 11 % 1000;
		new_points[i][1] = (n + i) * 17 % 1000;
		new_points[i].aug.id = n + i;
	});

	tree.BatchInsert(parlay::make_slice(new_points));
	std::cout << "Inserted " << insert_count << " new points" << std::endl;

	// 6. Batch delete some points
	size_t delete_count = 50;
	Points points_to_delete = points.subseq(0, delete_count);

	tree.BatchDelete(parlay::make_slice(points_to_delete));
	std::cout << "Deleted " << delete_count << " points" << std::endl;

	// 7. Clean up
	tree.DeleteTree();
	std::cout << "Example completed successfully!" << std::endl;
}

} // namespace kd_tree_example
