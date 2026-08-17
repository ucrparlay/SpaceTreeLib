#pragma once

#include <iostream>
#include <vector>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/p_tree.h"

// Example usage of p_tree from SpaceTreeLib
// p_tree uses space-filling curves (Morton/Hilbert) to order points

namespace p_tree_example
{

using coord_type = long;

// Define augmentation structure for points with space-filling curve code
// WARN: All functions must be defined
// p_tree requires aug_id_code which includes both an id and a curve code
struct aug_id_code {
	using id_type = int_fast32_t;
	using curve_code_type = uint64_t;

	aug_id_code() : code(0), id(0)
	{
	}

	void set_member(curve_code_type const &val)
	{
		code = val;
	}

	bool operator<(aug_id_code const &rhs) const
	{
		return code == rhs.code ? id < rhs.id : code < rhs.code;
	}

	bool operator==(aug_id_code const &rhs) const
	{
		// WARN: code is not important, we only need to ensure the id
		return id == rhs.id;
	}

	friend std::ostream &operator<<(std::ostream &os,
					aug_id_code const &rhs)
	{
		os << rhs.code << " " << rhs.id;
		return os;
	}

	curve_code_type code;
	id_type id;
};

// Define point type: 2D points with augmented ID and curve code
using Point = psi::aug_point<coord_type, 2, aug_id_code>;
using points_type = parlay::sequence<Point>;
using base_type = psi::point_traits<Point>;

// Define split rule using space-filling curve (Morton curve in this example)
using SplitRule = psi::spatial_filling_curve<psi::morton_curve<Point>>;

// Alternative: hilbert_curve<Point>
using another_split_rule_type =
	psi::spatial_filling_curve<psi::hilbert_curve<Point>>;

// Define p_tree type
// p_tree doesn't use LeafAug/InteriorAug like kd_tree/orth_tree
// It uses CPAM (Compressed Purely Functional Augmented Maps) internally to
// maintain the bounding box as an augmented value.
using Tree = psi::p_tree<psi::tree_traits<Point, another_split_rule_type>>;

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
		// Note: curve code will be computed during build
	});

	std::cout << "Created " << n << " random 2D points" << std::endl;

	// 2. build the p_tree
	Tree tree;
	tree.build(parlay::make_slice(points));
	std::cout << "Built PTree with " << n << " points using Morton curve"
		  << std::endl;

	// 3. K-Nearest neighbors query
	expect(tree.get_size() == n, "size after build");

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

	size_t brute = 0;
	for (auto const &p : points) {
		if (p[0] >= 400 && p[0] <= 600 && p[1] >= 400 && p[1] <= 600)
			brute++;
	}
	expect(count == brute, "range_query count matches brute force");
	expect(tree.range_count(query_box).first == brute,
	       "range_count matches brute force");

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

} // namespace p_tree_example
