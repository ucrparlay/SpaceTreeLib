/*
 * The code in docs/MANUAL.md, compiled and run. The manual went stale once
 * because nothing checked it; this is what stops that happening again. Keep
 * the snippets here identical to the ones in the document.
 */

#include <parlay/primitives.h>

#include <functional>
#include <ostream>

#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"
#include "psi/orth_tree.h"
#include "psi/p_tree.h"
#include "tests/unit/check.h"

/* --- "A point type" --- */
using coord_type = long;
using point = psi::basic_point<coord_type, 2>;

/* --- "A split rule" --- */
using split_rule = psi::orthogonal_split_rule<psi::rotate_dim<point>,
					      psi::spatial_median<point>>;

/* --- "Choosing and configuring a tree": the SPaC-tree point --- */
struct curve_code {
	using id_type = int_fast32_t;
	using curve_code_type = uint64_t;
	curve_code() : code(0), id(0)
	{
	}
	void set_member(curve_code_type const &v)
	{
		code = v;
	}
	bool operator<(curve_code const &r) const
	{
		return code == r.code ? id < r.id : code < r.code;
	}
	bool operator==(curve_code const &r) const
	{
		return id == r.id;
	}
	friend std::ostream &operator<<(std::ostream &o, curve_code const &r)
	{
		return o << r.code << " " << r.id;
	}
	curve_code_type code;
	id_type id;
};
using aug_pt = psi::aug_point<coord_type, 2, curve_code>;
using curve = psi::spatial_filling_curve<psi::morton_curve<aug_pt>>;

static parlay::sequence<point> sample_points(size_t n)
{
	parlay::sequence<point> points(n);
	for (size_t i = 0; i < n; i++) {
		points[i][0] = static_cast<coord_type>((i * 7) % 1000);
		points[i][1] = static_cast<coord_type>((i * 13) % 1000);
	}
	return points;
}

int main()
{
	CASE("manual/build");
	/* --- "Build" --- */
	psi::kd_tree<point, split_rule> tree;

	parlay::sequence<point> points = sample_points(1000);
	tree.build(parlay::make_slice(points));
	CHECK_EQ(tree.get_size(), points.size());

	CASE("manual/query-simple");
	/* --- "Query": the allocating forms --- */
	{
		typename decltype(tree)::box_type box;
		box.first[0] = 400;
		box.first[1] = 400;
		box.second[0] = 600;
		box.second[1] = 600;

		point query_point;
		query_point[0] = 500;
		query_point[1] = 500;

		parlay::sequence<point> inside = tree.range_query(box);
		auto nearest = tree.knn(query_point, 10);

		CHECK_EQ(nearest.size(), size_t(10));
		CHECK(nearest.front().second <= nearest.back().second);
		CHECK_EQ(inside.size(), tree.range_count(box).first);
		CHECK(!tree.empty());
	}

	CASE("manual/range");
	/* --- "Query": the buffer-taking forms --- */
	typename decltype(tree)::box_type box;
	box.first[0] = 400;
	box.first[1] = 400;
	box.second[0] = 600;
	box.second[1] = 600;

	auto [count, count_log] = tree.range_count(box);

	parlay::sequence<point> found(tree.get_size());
	auto [written, query_log] =
		tree.range_query(box, parlay::make_slice(found));
	CHECK_EQ(written, count);

	size_t brute = 0;
	for (auto const &p : points) {
		if (p[0] >= 400 && p[0] <= 600 && p[1] >= 400 && p[1] <= 600)
			brute++;
	}
	CHECK_EQ(count, brute);

	CASE("manual/knn");
	/* --- the knn snippet --- */
	using dis_type = typename point::dis_type;
	using nn_pair = std::pair<std::reference_wrapper<point>, dis_type>;

	point query_point;
	query_point[0] = 500;
	query_point[1] = 500;

	size_t const k = 10;
	parlay::sequence<nn_pair> result(k, nn_pair(std::ref(points[0]), 0));
	psi::bounded_queue<point, nn_pair> queue(parlay::make_slice(result));

	tree.knn(query_point, queue);
	CHECK_EQ(queue.size(), k);

	CASE("manual/update");
	/* --- "Update" --- */
	parlay::sequence<point> more_points = sample_points(200);
	for (auto &p : more_points)
		p[0] += 2000;
	auto points_to_remove = more_points;
	auto maybe_present = more_points;

	tree.batch_insert(parlay::make_slice(more_points));
	CHECK_EQ(tree.get_size(), points.size() + more_points.size());
	tree.batch_delete(parlay::make_slice(points_to_remove));
	CHECK_EQ(tree.get_size(), points.size());
	tree.batch_diff(parlay::make_slice(maybe_present));
	CHECK_EQ(tree.get_size(), points.size());

	CASE("manual/flatten");
	/* --- "Read everything back" --- */
	parlay::sequence<point> all(tree.get_size());
	size_t n = tree.flatten(parlay::make_slice(all));
	CHECK_EQ(n, tree.get_size());

	CASE("manual/parallel");
	/* --- "Parallelism" --- */
	auto queries = sample_points(64);
	parlay::sequence<size_t> counts(queries.size());
	parlay::parallel_for(0, queries.size(), [&](size_t i) {
		typename decltype(tree)::box_type b;
		for (size_t d = 0; d < point::get_dim(); d++) {
			b.first[d] = queries[i][d];
			b.second[d] = queries[i][d] + 50;
		}
		counts[i] = tree.range_count(b).first;
	});
	CHECK(counts.size() == queries.size());

	CASE("manual/trees");
	/* --- the three instantiations --- */
	psi::kd_tree<point, split_rule> kd;
	psi::orth_tree<point, split_rule> orth;
	psi::p_tree<aug_pt, curve> p;

	kd.build(parlay::make_slice(points));
	orth.build(parlay::make_slice(points));
	CHECK_EQ(kd.get_size(), points.size());
	CHECK_EQ(orth.get_size(), points.size());

	parlay::sequence<aug_pt> aug_points(500);
	for (size_t i = 0; i < aug_points.size(); i++) {
		aug_points[i][0] = static_cast<coord_type>((i * 7) % 1000);
		aug_points[i][1] = static_cast<coord_type>((i * 13) % 1000);
		aug_points[i].aug.id = static_cast<curve_code::id_type>(i);
	}
	p.build(parlay::make_slice(aug_points));
	CHECK_EQ(p.get_size(), aug_points.size());

	CASE("manual/cleanup");
	/* --- "Clean up": delete_tree may be called more than once --- */
	tree.delete_tree();
	tree.delete_tree();
	CHECK_EQ(tree.get_size(), size_t(0));

	return psi_test_result();
}
