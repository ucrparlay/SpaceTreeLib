/*
 * Batch updates and flatten. flatten is what rebuild is built on and nothing
 * exercised it: every flatten call in tests/ is commented out.
 */

#include <parlay/primitives.h>

#include "tests/unit/check.h"
#include "tests/unit/trees.h"

using namespace psi_test;

template <typename Kind>
static std::vector<typename Kind::point_type>
flatten_all(typename Kind::tree_type &tree)
{
	using point_type = typename Kind::point_type;
	parlay::sequence<point_type> out(tree.get_size());
	size_t n = tree.flatten(parlay::make_slice(out));
	CHECK_EQ(n, tree.get_size());
	return sorted(std::vector<point_type>(out.begin(), out.begin() + n));
}

/* build -> flatten returns exactly what went in */
template <typename Kind>
static void round_trip(char const *tag)
{
	CASE(tag);
	auto src = gen<Kind>(4000, 0, 500, seed_from_env());
	auto pts = make_points<Kind>(src);
	typename Kind::tree_type tree;
	tree.build(parlay::make_slice(pts));

	std::vector<typename Kind::point_type> ref(pts.begin(), pts.end());
	CHECK(flatten_all<Kind>(tree) == sorted(ref));

	/* a wrong-sized buffer must write nothing rather than overrun */
	parlay::sequence<typename Kind::point_type> small(7);
	CHECK_EQ(tree.flatten(parlay::make_slice(small)), size_t(0));
}

/* insert then delete the same points leaves the tree as it was */
template <typename Kind>
static void insert_delete_identity(char const *tag)
{
	CASE(tag);
	uint64_t seed = seed_from_env();
	auto base = gen<Kind>(3000, 0, 500, seed);
	auto extra = gen<Kind>(1500, 0, 500, seed + 1);
	/* ids must not collide: p_tree treats them as identity */
	for (size_t i = 0; i < extra.size(); i++)
		if constexpr (requires { extra[i].aug.id; })
			extra[i].aug.id =
				decltype(extra[i].aug.id)(base.size() + i);

	auto pts = make_points<Kind>(base);
	auto add = make_points<Kind>(extra);
	for (size_t i = 0; i < add.size(); i++)
		if constexpr (requires { add[i].aug.id; })
			add[i].aug.id =
				decltype(add[i].aug.id)(base.size() + i);

	typename Kind::tree_type tree;
	tree.build(parlay::make_slice(pts));
	auto before = flatten_all<Kind>(tree);

	tree.batch_insert(parlay::make_slice(add));
	CHECK_EQ(tree.get_size(), base.size() + extra.size());

	tree.batch_delete(parlay::make_slice(add));
	CHECK_EQ(tree.get_size(), base.size());
	CHECK(flatten_all<Kind>(tree) == before);
}

/* batch_diff tolerates points that were never inserted */
template <typename Kind>
static void diff_tolerates_absent(char const *tag)
{
	CASE(tag);
	uint64_t seed = seed_from_env();
	auto base = gen<Kind>(2000, 0, 500, seed + 2);
	auto absent = gen<Kind>(600, 900, 1200, seed + 3);
	for (size_t i = 0; i < absent.size(); i++)
		if constexpr (requires { absent[i].aug.id; })
			absent[i].aug.id = decltype(absent[i].aug.id)(9000 + i);

	auto pts = make_points<Kind>(base);
	auto gone = make_points<Kind>(absent);
	for (size_t i = 0; i < gone.size(); i++)
		if constexpr (requires { gone[i].aug.id; })
			gone[i].aug.id = decltype(gone[i].aug.id)(9000 + i);

	typename Kind::tree_type tree;
	tree.build(parlay::make_slice(pts));
	tree.batch_diff(parlay::make_slice(gone));
	CHECK_EQ(tree.get_size(), base.size());
}

/* enough churn to trip the imbalance rebuild, with answers still right */
template <typename Kind>
static void rebuild_keeps_answers(char const *tag)
{
	using point_type = typename Kind::point_type;
	CASE(tag);
	uint64_t seed = seed_from_env();
	auto src = gen<Kind>(6000, 0, 400, seed + 4);
	auto pts = make_points<Kind>(src);
	typename Kind::tree_type tree;
	tree.build(parlay::make_slice(pts));

	std::vector<point_type> live(pts.begin(), pts.end());
	for (int round = 0; round < 6; round++) {
		size_t chunk = live.size() / 3;
		parlay::sequence<point_type> del(chunk);
		for (size_t i = 0; i < chunk; i++)
			del[i] = live[i];
		tree.batch_delete(parlay::make_slice(del));
		live.erase(live.begin(), live.begin() + chunk);
		CHECK_EQ(tree.get_size(), live.size());

		auto fresh = gen<Kind>(chunk, 0, 400, seed + 100 + round);
		for (size_t i = 0; i < fresh.size(); i++)
			if constexpr (requires { fresh[i].aug.id; })
				fresh[i].aug.id = decltype(fresh[i].aug.id)(
					100000 + round * 10000 + i);
		auto ins = make_points<Kind>(fresh);
		for (size_t i = 0; i < ins.size(); i++)
			if constexpr (requires { ins[i].aug.id; })
				ins[i].aug.id = decltype(ins[i].aug.id)(
					100000 + round * 10000 + i);
		tree.batch_insert(parlay::make_slice(ins));
		live.insert(live.end(), ins.begin(), ins.end());
		CHECK_EQ(tree.get_size(), live.size());
	}
	CHECK(flatten_all<Kind>(tree) == sorted(live));
}

template <typename Kind>
static void run(char const *tag)
{
	round_trip<Kind>((std::string(tag) + "/flatten").c_str());
	insert_delete_identity<Kind>(
		(std::string(tag) + "/insert-delete").c_str());
	diff_tolerates_absent<Kind>((std::string(tag) + "/diff").c_str());
	rebuild_keeps_answers<Kind>((std::string(tag) + "/rebuild").c_str());
}

int main()
{
	run<kd_kind<long, 2>>("kd/2d");
	run<kd_kind<long, 3>>("kd/3d");
	run<orth_kind<long, 2>>("orth/2d");
	run<p_kind<long, 2>>("p/2d");
	run<kd_stretch_kind<long, 2>>("kd-stretch/2d");
	run<p_hilbert_kind<long, 2>>("p-hilbert/2d");
	return psi_test_result();
}
