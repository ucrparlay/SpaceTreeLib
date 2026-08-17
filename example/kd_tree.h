#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"

// Example usage of kd_tree from SpaceTreeLib
// This demonstrates basic operations: build, knn, range_query, batch_insert,
// and batch_delete

namespace kd_tree_example
{

// Type for each coordinate
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
using base_type = psi::point_traits<Point>;

// leaf_type augmentation: stores bounding box

// interior_type node augmentation: stores bounding box and parallel build flag

// Define split rule: max stretch dimension + object median
using SplitRule = psi::orthogonal_split_rule<psi::max_stretch_dim<Point>,
					     psi::object_median<Point>>;

// Alternative split rule: rotate dimension + spatial median
using another_split_rule_type =
	psi::orthogonal_split_rule<psi::rotate_dim<Point>,
				   psi::spatial_median<Point>>;

// Define kd_tree type
using Tree = psi::kd_tree<psi::tree_traits<Point, SplitRule>>;

/*
 * The example doubles as a smoke test: without these it prints success even
 * when every answer is wrong.
 */
static int example_failures = 0;

static void expect(bool ok, char const *what)
{
	if (!ok) {
		std::cerr << "  CHECK FAILED: " << what << std::endl;
		example_failures++;
	}
}

int run_example()
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

	// 2. build the kd_tree
	Tree tree;
	tree.build(parlay::make_slice(points));
	std::cout << "Built KdTree with " << n << " points" << std::endl;
	expect(tree.get_size() == n, "size after build");

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

	// 4. Range query (Range count is available as well)
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

	// the same answer by brute force, and range_count must agree
	size_t brute = 0;
	for (auto const &p : points) {
		if (p[0] >= 400 && p[0] <= 600 && p[1] >= 400 && p[1] <= 600)
			brute++;
	}
	expect(count == brute, "range_query count matches brute force");
	expect(tree.range_count(query_box).first == brute,
	       "range_count matches brute force");

	// the k-th nearest distance, by brute force
	std::vector<dis_type> all_dist;
	for (auto const &p : points) {
		dis_type dx = p[0] - query_point[0];
		dis_type dy = p[1] - query_point[1];
		all_dist.push_back(dx * dx + dy * dy);
	}
	std::sort(all_dist.begin(), all_dist.end());
	dis_type worst = 0;
	for (auto const &[pt, dist] : knn_result)
		worst = std::max(worst, dist);
	expect(worst == all_dist[K - 1],
	       "knn k-th distance matches brute force");

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
	expect(tree.get_size() == n + insert_count, "size after insert");

	// 6. Batch delete some points
	size_t delete_count = 50;
	points_type points_to_delete = points.subseq(0, delete_count);

	tree.batch_delete(parlay::make_slice(points_to_delete));
	std::cout << "Deleted " << delete_count << " points" << std::endl;
	expect(tree.get_size() == n + insert_count - delete_count,
	       "size after delete");

	// 7. Clean up
	tree.delete_tree();
	if (example_failures != 0) {
		std::cerr << "Example FAILED: " << example_failures
			  << " check(s)" << std::endl;
		return 1;
	}
	std::cout << "Example completed successfully!" << std::endl;
	return 0;
}

} // namespace kd_tree_example
