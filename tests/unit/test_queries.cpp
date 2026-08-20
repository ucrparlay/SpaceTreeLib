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
	return psi_test_result();
}
