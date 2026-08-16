/*
 * The allocating forms of knn and range_query, checked against the
 * buffer-taking ones they wrap.
 */

#include <parlay/primitives.h>

#include "tests/unit/check.h"
#include "tests/unit/trees.h"

using namespace psi_test;

template <typename Kind>
static void convenience_matches_primary(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using dis_type = typename point_type::dis_type;
	CASE(tag);

	uint64_t seed = seed_from_env();
	auto pts = make_points<Kind>(gen<Kind>(4000, 0, 300, seed));
	tree_type tree;
	CHECK(tree.empty());
	tree.build(parlay::make_slice(pts));
	CHECK(!tree.empty());

	std::vector<point_type> ref(pts.begin(), pts.end());

	/* range_query(box) must agree with range_query(box, buffer) */
	for (int t = 0; t < 8; t++) {
		auto corners = random_points<typename point_type::bp_type>(
			2, -20, 320, seed + 500 + t);
		typename tree_type::box_type box{corners[0], corners[1]};
		for (size_t d = 0; d < point_type::get_dim(); d++) {
			if (box.first.pnt[d] > box.second.pnt[d])
				std::swap(box.first.pnt[d], box.second.pnt[d]);
		}

		auto got = tree.range_query(box);
		CHECK_EQ(got.size(), tree.range_count(box).first);
		CHECK(sorted(std::vector<point_type>(got.begin(), got.end())) ==
		      range_query(ref, box));
	}

	/* knn(q, k) must agree with the queue form, and come back sorted */
	for (size_t k : {size_t(1), size_t(8), size_t(64)}) {
		auto queries = gen<Kind>(4, 0, 300, seed + 900);
		for (auto const &q : queries) {
			auto got = tree.knn(q, k);
			CHECK_EQ(got.size(), std::min(k, ref.size()));

			std::vector<dis_type> d;
			for (auto const &e : got)
				d.push_back(e.second);
			CHECK(std::is_sorted(d.begin(), d.end()));
			CHECK(near(d, knn_distances(ref, q, got.size())));
		}
	}

	/* degenerate: k of zero, and an empty tree */
	CHECK_EQ(tree.knn(ref[0], 0).size(), size_t(0));
	tree_type fresh;
	CHECK(fresh.empty());
	CHECK_EQ(fresh.knn(ref[0], 5).size(), size_t(0));
	typename tree_type::box_type any{ref[0], ref[0]};
	CHECK_EQ(fresh.range_query(any).size(), size_t(0));
}

int main()
{
	convenience_matches_primary<kd_kind<long, 2>>("kd/2d");
	convenience_matches_primary<kd_kind<long, 3>>("kd/3d");
	convenience_matches_primary<orth_kind<long, 2>>("orth/2d");
	convenience_matches_primary<p_kind<long, 2>>("p/2d");
	return psi_test_result();
}
