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

/* size, clear, bounds, count, contains, insert and erase. */
template <typename Kind>
static void container_api(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using box_type = typename tree_type::box_type;
	CASE(tag);

	uint64_t seed = seed_from_env();
	auto pts = make_points<Kind>(gen<Kind>(2000, 0, 200, seed));
	std::vector<point_type> ref(pts.begin(), pts.end());

	tree_type tree;
	CHECK(tree.empty());
	CHECK_EQ(tree.size(), size_t(0));
	tree.build(parlay::make_slice(pts));
	CHECK_EQ(tree.size(), ref.size());
	CHECK_EQ(tree.size(), tree.get_size());

	/* bounds covers every point, and is tight on a freshly built tree */
	box_type want{ref[0], ref[0]};
	for (auto const &p : ref) {
		for (size_t d = 0; d < point_type::get_dim(); d++) {
			want.first.pnt[d] = std::min(want.first.pnt[d], p[d]);
			want.second.pnt[d] = std::max(want.second.pnt[d], p[d]);
		}
	}
	box_type const &got = tree.bounds();
	for (size_t d = 0; d < point_type::get_dim(); d++) {
		CHECK_EQ(got.first.pnt[d], want.first.pnt[d]);
		CHECK_EQ(got.second.pnt[d], want.second.pnt[d]);
	}

	/* count is over coordinates, so brute force counts the same way */
	for (size_t i = 0; i < 12; i++) {
		auto const &p = ref[(i * 137) % ref.size()];
		CHECK_EQ(tree.count(p), range_count(ref, box_type(p, p)));
		CHECK(tree.contains(p));
	}

	/*
	 * A point inside bounds that is not there yet -- orth_tree refuses to
	 * grow its box, so an outside point would be a precondition failure.
	 */
	point_type fresh_point = ref[0];
	bool found_gap = false;
	for (typename point_type::coord_type c = 0; c < 200; c++) {
		for (size_t d = 0; d < point_type::get_dim(); d++)
			fresh_point.pnt[d] = want.first.pnt[d] + c;
		if constexpr (requires { fresh_point.aug.id; })
			fresh_point.aug.id = 1 << 20;
		if (!tree.contains(fresh_point)) {
			found_gap = true;
			break;
		}
	}
	CHECK(found_gap);

	size_t const n = tree.size();
	tree.insert(fresh_point);
	CHECK_EQ(tree.size(), n + 1);
	CHECK(tree.contains(fresh_point));

	CHECK_EQ(tree.erase(fresh_point), size_t(1));
	CHECK_EQ(tree.size(), n);
	CHECK(!tree.contains(fresh_point));
	/* erasing what is not there removes nothing and says so */
	CHECK_EQ(tree.erase(fresh_point), size_t(0));
	CHECK_EQ(tree.size(), n);

	tree.clear();
	CHECK(tree.empty());
	CHECK_EQ(tree.size(), size_t(0));
	CHECK(!tree.contains(ref[0]));
	CHECK_EQ(tree.count(ref[0]), size_t(0));
}

int main()
{
	convenience_matches_primary<kd_kind<long, 2>>("kd/2d");
	convenience_matches_primary<kd_kind<long, 3>>("kd/3d");
	convenience_matches_primary<orth_kind<long, 2>>("orth/2d");
	convenience_matches_primary<p_kind<long, 2>>("p/2d");

	container_api<kd_kind<long, 2>>("kd/2d/container");
	container_api<kd_kind<long, 3>>("kd/3d/container");
	container_api<orth_kind<long, 2>>("orth/2d/container");
	container_api<p_kind<long, 2>>("p/2d/container");
	return psi_test_result();
}
