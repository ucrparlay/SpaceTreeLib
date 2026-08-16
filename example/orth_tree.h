#pragma once

#include <iostream>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/orth_tree.h"

// Example usage of orth_tree from SpaceTreeLib
// orth_tree is an orthogonal tree that splits space into 2^d regions at each
// node (quadtree for 2D, octree for 3D)

namespace orth_tree_example
{

using coord_type = long;

// Define augmentation structure for points (stores an ID)
// WARN: All functions must be defined
struct aug_id {
	using id_type = int;
	id_type id;

	bool operator<(aug_id const &rhs) const
	{
		return id < rhs.id;
	}
	bool operator==(aug_id const &rhs) const
	{
		return id == rhs.id;
	}
	friend std::ostream &operator<<(std::ostream &os, aug_id const &rhs)
	{
		os << rhs.id;
		return os;
	}
};

// Define point type: 2D points with augmented ID
using Point = psi::aug_point<coord_type, 2, aug_id>;
using points_type = parlay::sequence<Point>;
using base_type = psi::base_tree<Point>;

// leaf_type augmentation: stores nothing

// leaf_type augmentation: stores bounding box

// interior_type node augmentation: stores bounding box and parallel build flag
// interior_type node augmentation with box

// Helper function to determine skeleton height for orth_tree
// Same as in test_framework.h
// Define split rule: stadard rotate dimension + spatial median
using SplitRule = psi::orthogonal_split_rule<psi::rotate_dim<Point>,
					     psi::spatial_median<Point>>;

// Define orth_tree type with augmentations
using Tree = psi::orth_tree<Point, SplitRule, psi::box_leaf_aug<base_type>,
			    psi::box_interior_aug<base_type>>;

void run_example()
{
	// 1. create sample 2D points
	size_t n = 1000;
	points_type points(n);

	// Generate random points in a 1000x1000 grid
	parlay::parallel_for(0, n, [&](size_t i) {
		points[i][0] = (i * 7) % 1000;	// x coordinate
		points[i][1] = (i * 13) % 1000; // y coordinate
		points[i].aug.id = i;		// unique ID
	});

	std::cout << "Created " << n << " random 2D points" << std::endl;

	// 2. build the orth_tree
	// orth_tree requires a bounding box to be specified
	Tree tree;

	// Calculate bounding box for the points
	// NOTE: This is not necessary, unless you want to specify a certain
	// bounding box for the tree.
	typename Tree::box_type bounding_box =
		Tree::get_box(parlay::make_slice(points));

	std::cout << "Bounding box: [(" << bounding_box.first[0] << ", "
		  << bounding_box.first[1] << "), (" << bounding_box.second[0]
		  << ", " << bounding_box.second[1] << ")]" << std::endl;

	tree.build(parlay::make_slice(points), bounding_box);
	std::cout << "Built OrthTree (quadtree) with " << n << " points"
		  << std::endl;

	// 3. K-Nearest neighbors query
	int K = 10;
	Point query_point;
	query_point[0] = 500;
	query_point[1] = 500;
	query_point.aug.id = -1;

	using dis_type = typename Point::dis_type;
	using nn_pair = std::pair<std::reference_wrapper<Point>, dis_type>;

	parlay::sequence<nn_pair> knn_result(K,
					     nn_pair(std::ref(points[0]), 0));
	psi::bounded_queue<Point, nn_pair> bq(parlay::make_slice(knn_result));

	tree.knn(query_point, bq);

	std::cout << "Found " << K
		  << " nearest neighbors to point (500, 500) (unsorted):"
		  << std::endl;
	for (int i = 0; i < std::min(3, K); i++) {
		auto &[pt, dist] = knn_result[i];
		std::cout << "  Point " << pt.get().aug.id << " at ("
			  << pt.get()[0] << ", " << pt.get()[1]
			  << ") with distance " << dist << std::endl;
	}

	// 4. Range query
	typename Tree::box_type query_box;
	query_box.first[0] = 400;
	query_box.first[1] = 400;
	query_box.second[0] = 600;
	query_box.second[1] = 600;

	points_type range_result(n); // Allocate max possible size
	auto [count, logger] =
		tree.range_query(query_box, parlay::make_slice(range_result));

	std::cout << "Range query [400,600]x[400,600] found " << count
		  << " points" << std::endl;

	// 5. Batch insert new points
	size_t insert_count = 100;
	points_type new_points(insert_count);
	parlay::parallel_for(0, insert_count, [&](size_t i) {
		new_points[i][0] = (n + i) * 11 % 1000;
		new_points[i][1] = (n + i) * 17 % 1000;
		new_points[i].aug.id = n + i;
	});

	tree.batch_insert(parlay::make_slice(new_points));
	std::cout << "Inserted " << insert_count << " new points" << std::endl;

	// 6. Batch delete some points
	size_t delete_count = 50;
	points_type points_to_delete = points.subseq(0, delete_count);

	tree.batch_delete(parlay::make_slice(points_to_delete));
	std::cout << "Deleted " << delete_count << " points" << std::endl;

	// 7. Clean up
	tree.delete_tree();
	std::cout << "Example completed successfully!" << std::endl;
}

} // namespace orth_tree_example
