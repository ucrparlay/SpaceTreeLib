/*
 * The inputs a user reaches by accident: an unbuilt tree, k outside the tree's
 * size, a query box that covers nothing, an output buffer that is too small,
 * and floating-point coordinates.
 */

#include <parlay/primitives.h>

#include "tests/unit/check.h"
#include "tests/unit/trees.h"

using namespace psi_test;

template <typename Kind>
static void unbuilt_tree(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	CASE(tag);

	tree_type tree;
	CHECK_EQ(tree.get_size(), size_t(0));

	typename point_type::coord_type lo[point_type::get_dim()];
	typename point_type::coord_type hi[point_type::get_dim()];
	for (size_t d = 0; d < point_type::get_dim(); d++) {
		lo[d] = 0;
		hi[d] = 100;
	}
	typename tree_type::box_type box{typename point_type::bp_type(lo),
					 typename point_type::bp_type(hi)};

	CHECK_EQ(tree.range_count(box).first, size_t(0));
	parlay::sequence<point_type> out(0);
	CHECK_EQ(tree.range_query(box, parlay::make_slice(out)).first,
		 size_t(0));
	CHECK_EQ(tree.flatten(parlay::make_slice(out)), size_t(0));
}

template <typename Kind>
static void empty_and_inverted_box(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	CASE(tag);

	auto pts = make_points<Kind>(gen<Kind>(2000, 0, 200, seed_from_env()));
	tree_type tree;
	tree.build(parlay::make_slice(pts));

	/* inverted in one dimension: covers nothing */
	typename point_type::coord_type lo[point_type::get_dim()];
	typename point_type::coord_type hi[point_type::get_dim()];
	for (size_t d = 0; d < point_type::get_dim(); d++) {
		lo[d] = 0;
		hi[d] = 200;
	}
	lo[0] = 150;
	hi[0] = 50;
	typename tree_type::box_type inverted{typename point_type::bp_type(lo),
					      typename point_type::bp_type(hi)};
	CHECK_EQ(tree.range_count(inverted).first, size_t(0));
	parlay::sequence<point_type> out(pts.size());
	CHECK_EQ(tree.range_query(inverted, parlay::make_slice(out)).first,
		 size_t(0));

	/* a single point is a legal, non-empty box */
	for (size_t d = 0; d < point_type::get_dim(); d++) {
		lo[d] = 7;
		hi[d] = 7;
	}
	typename tree_type::box_type dot{typename point_type::bp_type(lo),
					 typename point_type::bp_type(hi)};
	std::vector<point_type> ref(pts.begin(), pts.end());
	CHECK_EQ(tree.range_count(dot).first, range_count(ref, dot));
}

template <typename Kind>
static void undersized_output(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	CASE(tag);

	auto pts = make_points<Kind>(gen<Kind>(3000, 0, 100, seed_from_env()));
	tree_type tree;
	tree.build(parlay::make_slice(pts));

	typename point_type::coord_type lo[point_type::get_dim()];
	typename point_type::coord_type hi[point_type::get_dim()];
	for (size_t d = 0; d < point_type::get_dim(); d++) {
		lo[d] = 0;
		hi[d] = 100;
	}
	typename tree_type::box_type all{typename point_type::bp_type(lo),
					 typename point_type::bp_type(hi)};
	size_t total = tree.range_count(all).first;
	CHECK(total > 0);

	for (size_t cap : {size_t(0), size_t(1), size_t(50), total / 2}) {
		parlay::sequence<point_type> out(cap);
		size_t got =
			tree.range_query(all, parlay::make_slice(out)).first;
		CHECK_EQ(got, cap);
	}
}

template <typename Kind>
static void k_outside_range(char const *tag)
{
	using point_type = typename Kind::point_type;
	using tree_type = typename Kind::tree_type;
	using dis_type = typename point_type::dis_type;
	CASE(tag);

	auto pts = make_points<Kind>(gen<Kind>(500, 0, 100, seed_from_env()));
	tree_type tree;
	tree.build(parlay::make_slice(pts));
	std::vector<point_type> ref(pts.begin(), pts.end());
	auto q = ref[3];

	using nn_pair = std::pair<std::reference_wrapper<point_type>, dis_type>;

	/* k == 0: nothing to fill, and nothing outside the range written */
	{
		parlay::sequence<nn_pair> buf(0, nn_pair(std::ref(pts[0]), 0));
		psi::bounded_queue<point_type, nn_pair> bq(
			parlay::make_slice(buf));
		tree.knn(q, bq);
		CHECK_EQ(bq.size(), size_t(0));
	}
	/* k larger than the tree */
	{
		size_t k = ref.size() + 25;
		parlay::sequence<nn_pair> buf(k, nn_pair(std::ref(pts[0]), 0));
		psi::bounded_queue<point_type, nn_pair> bq(
			parlay::make_slice(buf));
		tree.knn(q, bq);
		/* only the filled prefix is meaningful; the rest of the queue
		 * still holds whatever the caller put there */
		CHECK_EQ(bq.size(), ref.size());
		std::vector<dis_type> got;
		for (size_t i = 0; i < bq.size(); i++)
			got.push_back(buf[i].second);
		std::sort(got.begin(), got.end());
		CHECK(near(got, knn_distances(ref, q, ref.size())));
	}
}

template <typename Kind>
static void run(char const *tag)
{
	unbuilt_tree<Kind>((std::string(tag) + "/unbuilt").c_str());
	empty_and_inverted_box<Kind>((std::string(tag) + "/empty-box").c_str());
	undersized_output<Kind>((std::string(tag) + "/short-buffer").c_str());
	k_outside_range<Kind>((std::string(tag) + "/k-range").c_str());
}

int main()
{
	run<kd_kind<long, 2>>("kd/2d");
	run<orth_kind<long, 2>>("orth/2d");
	run<p_kind<long, 2>>("p/2d");
	/* double coordinates compile a comparison branch that the integer
	 * instantiations never reach */
	run<kd_kind<double, 2>>("kd/2d-double");
	return psi_test_result();
}
