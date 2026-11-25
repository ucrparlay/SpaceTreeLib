/*
 * Queries against brute force: range_count, range_query and knn on every
 * tree, over ordinary and degenerate inputs.
 */

#include <parlay/primitives.h>

#include "tests/unit/check.h"
#include "tests/unit/trees.h"

using namespace psi_test;

template <typename Kind>
static void check_queries(std::vector<typename Kind::point_type> const &src,
			  char const *what)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using box_type = typename tree_type::box_type;

	CASE(what);
	auto pts = make_points<Kind>(src);
	tree_type tree;
	tree.build(parlay::make_slice(pts));
	CHECK_EQ(tree.get_size(), src.size());

	std::vector<point_type> ref(pts.begin(), pts.end());
	uint64_t seed = seed_from_env();

	/* boxes drawn from the data, so they are neither always empty nor
	 * always everything */
	for (int t = 0; t < 12; t++) {
		auto corners = random_points<typename point_type::bp_type>(
			2, -40, 140, seed + 1000 + t);
		box_type box{corners[0], corners[1]};
		for (size_t d = 0; d < point_type::get_dim(); d++) {
			if (box.first.pnt[d] > box.second.pnt[d])
				std::swap(box.first.pnt[d], box.second.pnt[d]);
		}
		size_t expect = range_count(ref, box);
		CHECK_EQ(tree.range_count(box).first, expect);

		parlay::sequence<point_type> out(ref.size());
		size_t got =
			tree.range_query(box, parlay::make_slice(out)).first;
		CHECK_EQ(got, expect);

		std::vector<point_type> got_pts(out.begin(), out.begin() + got);
		CHECK(sorted(got_pts) == range_query(ref, box));
	}

	/* knn: compare distances, since ties leave the k-th neighbour's
	 * identity ambiguous but not its distance */
	using dis_type = typename point_type::dis_type;
	for (size_t k : {size_t(1), size_t(5), size_t(32), ref.size() + 7}) {
		if (ref.empty())
			break;
		auto queries = gen<Kind>(6, -20, 120, seed + 77);
		for (auto const &q : queries) {
			size_t kk = std::min(k, ref.size());
			using nn_pair =
				std::pair<std::reference_wrapper<point_type>,
					  dis_type>;
			parlay::sequence<nn_pair> buf(
				kk, nn_pair(std::ref(pts[0]), 0));
			psi::bounded_queue<point_type, nn_pair> bq(
				parlay::make_slice(buf));
			tree.knn(q, bq);

			std::vector<dis_type> got;
			for (size_t i = 0; i < bq.size(); i++)
				got.push_back(buf[i].second);
			std::sort(got.begin(), got.end());
			CHECK(near(got, knn_distances(ref, q, kk)));
		}
	}
}

template <typename Kind>
static void run(char const *tag)
{
	uint64_t seed = seed_from_env();
	check_queries<Kind>(gen<Kind>(3000, 0, 100, seed),
			    (std::string(tag) + "/random").c_str());
	check_queries<Kind>(gen<Kind>(1, 0, 100, seed),
			    (std::string(tag) + "/one-point").c_str());
	check_queries<Kind>(
		[] {
			auto v = repeated_points<typename Kind::point_type>(500,
									    42);
			for (size_t i = 0; i < v.size(); i++)
				if constexpr (requires { v[i].aug.id; })
					v[i].aug.id = decltype(v[i].aug.id)(i);
			return v;
		}(),
		(std::string(tag) + "/all-identical").c_str());
	/* every point on one plane: the split rules have to cope with a box
	 * that is degenerate in a dimension */
	check_queries<Kind>(
		[&] {
			auto v = gen<Kind>(800, 0, 100, seed + 5);
			for (auto &p : v)
				p.pnt[0] = 50;
			return v;
		}(),
		(std::string(tag) + "/one-plane").c_str());
}

/* Two clusters at opposite corners of [0, 2^31)^3: the squared distance
 * between them passes INT64_MAX. */
template <typename Kind>
static void far_clusters(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using coord_type = typename point_type::coord_type;
	using dis_type = typename point_type::dis_type;
	using nn_pair = std::pair<std::reference_wrapper<point_type>, dis_type>;
	CASE(tag);

	uint64_t seed = seed_from_env();
	coord_type const top = (coord_type(1) << 31) - 1;
	auto low = gen<Kind>(500, 0, 100, seed);
	auto high = gen<Kind>(500, top - 100, top, seed + 1);
	std::vector<point_type> src(low);
	src.insert(src.end(), high.begin(), high.end());
	auto pts = make_points<Kind>(src);
	tree_type tree;
	tree.build(parlay::make_slice(pts));

	/* the neighbours of a point are all in its own cluster */
	size_t const k = 10;
	for (size_t i = 0; i < src.size(); i += 25) {
		auto const &own = i < low.size() ? low : high;
		parlay::sequence<nn_pair> buf(k, nn_pair(std::ref(pts[0]), 0));
		psi::bounded_queue<point_type, nn_pair> bq(
			parlay::make_slice(buf));
		tree.knn(src[i], bq);

		std::vector<dis_type> got;
		for (size_t j = 0; j < bq.size(); j++)
			got.push_back(buf[j].second);
		std::sort(got.begin(), got.end());
		CHECK(near(got, knn_distances(own, src[i], k)));
	}
}

/*
 * Coordinates past the bits a code holds an axis. With the box, the code
 * holds their top bits, so a leaf holds near points and a knn reads few
 * leaves.
 */
template <typename Kind>
static void top_bits(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using box_type = typename tree_type::box_type;
	using bp_type = typename point_type::bp_type;
	using coord_type = typename point_type::coord_type;
	using dis_type = typename point_type::dis_type;
	using nn_pair = std::pair<std::reference_wrapper<point_type>, dis_type>;
	CASE(tag);

	uint64_t seed = seed_from_env();
	coord_type const top = (coord_type(1) << 31) - 1;
	auto src = gen<Kind>(20000, 0, top, seed);
	auto pts = make_points<Kind>(src);
	tree_type tree;
	tree.set_bounding_box(box_type(bp_type(coord_type(0)), bp_type(top)));
	tree.build(parlay::make_slice(pts));

	/* Another tree's update sets the shared shift to its own. */
	auto few = make_points<Kind>(gen<Kind>(100, 0, 100, seed + 1));
	tree_type other;
	other.build(parlay::make_slice(few));

	size_t const half = src.size() / 2;
	tree.batch_delete(parlay::make_slice(pts).cut(0, half));
	CHECK_EQ(tree.get_size(), src.size() - half);

	std::vector<point_type> const ref(src.begin() + half, src.end());
	size_t const k = 10;
	size_t leaf_cnt = 0, query_cnt = 0;
	for (size_t i = half; i < src.size(); i += 50) {
		parlay::sequence<nn_pair> buf(k, nn_pair(std::ref(pts[0]), 0));
		psi::bounded_queue<point_type, nn_pair> bq(
			parlay::make_slice(buf));
		leaf_cnt += tree.knn(src[i], bq).vis_leaf_num;
		query_cnt++;

		std::vector<dis_type> got;
		for (size_t j = 0; j < bq.size(); j++)
			got.push_back(buf[j].second);
		std::sort(got.begin(), got.end());
		CHECK(near(got, knn_distances(ref, src[i], k)));
	}
	CHECK(leaf_cnt < 30 * query_cnt);
}

/*
 * One far value on the last axis stretches the tree's box, so the cuts of
 * that axis leave the other points together and a leaf's part of the box is
 * far taller than its points. A box reads no leaf whose points lie apart from
 * it or all inside it.
 */
template <typename Kind>
static void far_value(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using box_type = typename tree_type::box_type;
	using bp_type = typename point_type::bp_type;
	using coord_type = typename point_type::coord_type;
	CASE(tag);

	uint64_t seed = seed_from_env();
	size_t const last = point_type::get_dim() - 1;
	auto src = gen<Kind>(20000, 0, 1000, seed);
	src.push_back(src[0]);
	src.back().pnt[last] = coord_type(1) << 30;
	auto pts = make_points<Kind>(src);
	tree_type tree;
	tree.build(parlay::make_slice(pts));
	std::vector<point_type> ref(pts.begin(), pts.end());

	/* checks the answers of a box over [low, high] on the last axis and
	 * returns the leaves it read */
	auto query = [&](coord_type low, coord_type high) {
		box_type box(bp_type(coord_type(0)), bp_type(coord_type(1000)));
		box.first.pnt[last] = low;
		box.second.pnt[last] = high;
		size_t const expect = range_count(ref, box);

		auto const counted = tree.range_count(box);
		CHECK_EQ(counted.first, expect);

		parlay::sequence<point_type> out(ref.size());
		auto const reported =
			tree.range_query(box, parlay::make_slice(out));
		CHECK_EQ(reported.first, expect);
		std::vector<point_type> got(out.begin(),
					    out.begin() + reported.first);
		CHECK(sorted(got) == range_query(ref, box));
		return counted.second.vis_leaf_num +
		       reported.second.vis_leaf_num;
	};
	auto check = [&] {
		/* between the near points and the far one */
		CHECK_EQ(query(5000, 6000), size_t(0));
		/* around the near points */
		CHECK_EQ(query(0, 2000), size_t(0));
		/* through them */
		for (coord_type low = 0; low < 1000; low += 97)
			query(low, low + 10);
	};
	check();

	size_t const half = 10000;
	tree.batch_delete(parlay::make_slice(pts).cut(0, half));
	ref.erase(ref.begin(), ref.begin() + half);
	check();

	auto more = make_points<Kind>(gen<Kind>(5000, 0, 1000, seed + 1));
	tree.batch_insert(parlay::make_slice(more));
	ref.insert(ref.end(), more.begin(), more.end());
	check();
}

int main()
{
	run<kd_kind<long, 2>>("kd/2d");
	run<kd_kind<long, 3>>("kd/3d");
	run<orth_kind<long, 2>>("orth/2d");
	run<orth_kind<long, 3>>("orth/3d");
	run<p_kind<long, 2>>("p/2d");
	run<kd_stretch_kind<long, 2>>("kd-stretch/2d");
	run<kd_stretch_kind<long, 3>>("kd-stretch/3d");
	run<p_hilbert_kind<long, 2>>("p-hilbert/2d");
	far_clusters<kd_kind<long, 3>>("kd/3d/far-clusters");
	far_clusters<orth_kind<long, 3>>("orth/3d/far-clusters");
	far_clusters<p_kind<long, 3>>("p/3d/far-clusters");
	far_clusters<kd_stretch_kind<long, 3>>("kd-stretch/3d/far-clusters");
	far_clusters<p_hilbert_kind<long, 3>>("p-hilbert/3d/far-clusters");
	top_bits<p_kind<long, 3>>("p/3d/top-bits");
	top_bits<p_hilbert_kind<long, 3>>("p-hilbert/3d/top-bits");
	far_value<orth_kind<long, 2>>("orth/2d/far-value");
	far_value<orth_kind<long, 3>>("orth/3d/far-value");
	return psi_test_result();
}
