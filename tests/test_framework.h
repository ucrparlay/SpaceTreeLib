#pragma once

#include <algorithm>
#include <filesystem>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <iostream>

#include "seed.h"
#include "baselines/cpam_raw/cpamtree.hpp"
#include "baselines/zdtree/zdtree.hpp"
#include "baselines/zdtree_3d/zdtree_3d.hpp"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "common/time_loop.h"
#include "psi/dependence/concepts.h"
#include "parlay/internal/group_by.h"
#include "parlay/monoid.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/slice.h"
#include "psi/base_tree.h"
#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"
#include "psi/orth_tree.h"
#include "psi/p_tree.h"

#ifdef CCP
using coord_type = long;
// using coord_type = double;
#else
using coord_type = long;
// using coord_type = double;
#endif // CCP

using scalar_type = coord_type;
using namespace psi;

// NOTE: knn size
static constexpr double batch_query_ratio = 0.01;
static constexpr size_t batch_query_osm_size = 10000000;
// NOTE: rectangle numbers
static constexpr int range_query_num = 50000;
static constexpr int single_query_log_repeat_num = 100;

// NOTE: rectangle numbers for inba ratio
static constexpr int rangeQueryNumInbaRatio = 50000;
// NOTE: insert batch ratio for inba ratio
static constexpr double insertBatchInbaRatio = 0.001;
// NOTE: knn batch ratio for inba ratio
static constexpr double knnBatchInbaRatio = 0.1;

// NOTE: Insert Ratio when summary
static constexpr double batch_insert_ratio = 0.01;
// NOTE: DIff Ratio when summary
static constexpr double batch_diff_total_ratio = 0.01;
static constexpr double batch_diff_overlap_ratio = 0.2;

// NOTE: rectange type used in summary
static constexpr int summaryRangeQueryType = 2;
// NOTE: range query num in summary
static constexpr int summary_range_query_num = 50000;

// NOTE: helper for delete type
enum delete_type { kBatchDelete, kBatchDiff };

// * [a,b)
inline size_t get_random_index(size_t a, size_t b, [[maybe_unused]] int seed)
{
	return size_t((rand() % (b - a)) + a);
}

template <typename Point, typename Tree, bool SavePoint>
size_t recurse_box(parlay::slice<Point *, Point *> in, auto &box_seq, int DIM,
		   std::pair<size_t, size_t> range, int &idx, int rec_num,
		   int type)
{
	size_t n = in.size();
	if (idx >= rec_num || n < range.first || n == 0)
		return 0;

	size_t mx = 0;
	bool goon = false;
	if (n <= range.second) {
		if constexpr (SavePoint) {
			box_seq[idx++] = std::make_pair(
				Tree::get_box(in), parlay::to_sequence(in));
		} else {
			box_seq[idx++] =
				std::make_pair(Tree::get_box(in), in.size());
		}

		// WARN: Modify the coefficient to make the rectangle size
		// distribute as uniform as possible within the range in
		// lograthmic scale
		if (parlay::all_of(
			    // if ((type == 0 && n < 40 * range.first) ||
			    // (type == 1 && n < 10 * range.first) ||
			    // (type == 2 && n < 2 * range.first) ||
			    // parlay::all_of(
			    in, [&](Point const &p) {
				    return p.same_dimension(in[0]);
			    })) {
			// NOTE: handle the cose that all points_type are the
			// same which is undivideable
			return in.size();
		} else {
			goon = true;
			mx = n;
		}
	}

	int dim = get_random_index(0, DIM, rand());
	size_t pos = get_random_index(0, n, rand());
	parlay::sequence<bool> flag(n, 0);
	parlay::parallel_for(0, n, [&](size_t i) {
		if (psi::num_comparator<coord_type>::gt(in[i][dim],
							in[pos][dim]))
			flag[i] = 1;
		else
			flag[i] = 0;
	});
	auto [out, m] = parlay::internal::split_two(in, flag);

	assert(out.size() == n);
	// std::cout << dim << " " << out[0] << out[m] << std::endl;
	size_t l = 0, r = 0;
	l = recurse_box<Point, Tree, SavePoint>(out.cut(0, m), box_seq, DIM,
						range, idx, rec_num, type);
	r = recurse_box<Point, Tree, SavePoint>(out.cut(m, n), box_seq, DIM,
						range, idx, rec_num, type);

	if (goon) {
		return mx;
	} else {
		return std::max(l, r);
	}
}

template <typename Point, typename Tree, bool SavePoint, bool FixSize = false>
auto gen_rectangles(int rec_num, int const type,
		    parlay::sequence<Point> const &WP, int DIM)
{
	using points_type = typename Tree::points_type;
	using box_type = typename Tree::box_type;
	using box_seq_type = std::conditional_t<
		SavePoint == false,
		parlay::sequence<std::pair<box_type, size_t>>,
		parlay::sequence<std::pair<box_type, points_type>>>;

	size_t n = WP.size();
	std::pair<size_t, size_t> range;
	if constexpr (FixSize) {
		if (type == 0) { //* small bracket
			range.first = 1;
			// range.second =
			// static_cast<size_t>(std::sqrt(std::sqrt(1.0 * n)));
			range.second = 100;
		} else if (type == 1) { //* medium bracket
			range.first = 100;
			// range.second = static_cast<size_t>(std::sqrt(1.0 *
			// n));
			range.second = 10000;
		} else if (type == 2) { //* large bracket
			range.first = 10000;
			if (n > 1'000'000) {
				range.second = 1'000'000;
			} // ensure we can generate 50k rect.
			else {
				range.second = n - 1;
			}
		}
	} else {
		if (type == 0) { //* small bracket
			range.first = 1;
			range.second = static_cast<size_t>(
				std::sqrt(std::sqrt(1.0 * n)));
		} else if (type == 1) { //* medium bracket
			range.second = static_cast<size_t>(
				std::sqrt(std::sqrt(1.0 * n)));
			range.second = static_cast<size_t>(std::sqrt(1.0 * n));
		} else if (type == 2) { //* large bracket
			range.second = static_cast<size_t>(std::sqrt(1.0 * n));
			if (n > 1'000'000) {
				range.second = 1'000'000;
			} // ensure we can generate 50k rect.
			else {
				range.second = n - 1;
			}
		}
	}
	/*
	 * Clamp to the input: the large bracket asks for 10,000-point
	 * rectangles, so below that recurse_box records nothing and the loop
	 * below never terminates. For the benchmark sizes the clamp never
	 * fires.
	 */
	range.first = std::min(range.first, n);
	range.second = std::max(range.first, std::min(range.second, n));

	box_seq_type box_seq(rec_num);
	int cnt = 0;
	points_type wp(n);

	srand(10);

	// std::cout << " " << range.first << " " << range.second << std::endl;

	size_t max_size = 0;
	int barren_passes = 0;
	while (cnt < rec_num) {
		int const cnt_before = cnt;
		parlay::copy(WP, wp);
		auto r = recurse_box<Point, Tree, SavePoint>(
			parlay::make_slice(wp), box_seq, DIM, range, cnt,
			rec_num, type);
		max_size = std::max(max_size, r);
		/* Bail out rather than spin forever if the input cannot
		 * produce the requested rectangles. */
		barren_passes = (cnt == cnt_before) ? barren_passes + 1 : 0;
		if (barren_passes >= 100) {
			std::cout << "gen_rectangles: only " << cnt << " of "
				  << rec_num << " rectangles fit in " << n
				  << " points; continuing with those\n";
			box_seq.resize(cnt > 0 ? cnt : 0);
			break;
		}
	}
	// std::cout << "finish generate " << std::endl;
	return std::make_pair(box_seq, max_size);
}

template <typename Point, typename Tree, bool test_time = true, int print = 1>
void build_tree(parlay::sequence<Point> const &WP, int const &rounds, Tree &pkd,
		int remaining_frac = 1)
{
	using points_type = typename Tree::points_type;
	using leaf_type = typename Tree::leaf_type;
	using interior_type = typename Tree::interior_type;

	double loop_late = rounds > 1 ? 0.01 : -100;
	size_t n = WP.size();
	// size_t n = 100;
	points_type wp = points_type::uninitialized(n);

	if constexpr (test_time) {
		pkd.delete_tree();

		double aveBuild = time_loop(
			rounds, loop_late,
			[&]() { parlay::copy(WP.cut(0, n), wp.cut(0, n)); },
			[&]() { pkd.build(wp.cut(0, n)); },
			[&]() { pkd.delete_tree(); });

		parlay::copy(WP.cut(0, n / remaining_frac),
			     wp.cut(0, n / remaining_frac));
		pkd.build(wp.cut(0, n / remaining_frac));

		if constexpr (print == 0) {
			std::cout << aveBuild << " " << std::flush;
		} else if constexpr (print == 1) {
			std::cout << aveBuild << " " << std::flush;
			if constexpr (is_kd_tree<Tree> || is_orth_tree<Tree>) {
				auto deep = pkd.template get_ave_tree_height<
					leaf_type, interior_type>();
				std::cout << deep << " " << std::flush;
			} else {
				std::cout << "-1"
					  << " " << std::flush;
			}
		} else if constexpr (print == 2) {
			size_t max_deep = 0;
			std::cout << aveBuild << " ";
			if constexpr (is_kd_tree<Tree> || is_orth_tree<Tree>) {
				std::cout << pkd.template get_max_tree_depth<
						     leaf_type, interior_type>(
						     pkd.get_root(), max_deep)
					  << " "
					  << pkd.template get_ave_tree_height<
						     leaf_type, interior_type>()
					  << " " << std::flush;
			} else {
				std::cout << "-1 -1"
					  << " " << std::flush;
			}
		} else if constexpr (print == 3) { // for incre insert directly
			puts("# Insert");
			std::cout << "## " << 1 << std::endl;
			std::cout << "median: (1, " << aveBuild
				  << ")-> min: (1, " << aveBuild
				  << ")-> max: (1, " << aveBuild
				  << ")-> tot: " << aveBuild
				  << "-> avg: " << aveBuild << std::endl;
		}

	} else {
		// NOTE: always return a built tree
		pkd.delete_tree();
		parlay::copy(WP.cut(0, n / remaining_frac),
			     wp.cut(0, n / remaining_frac));
		pkd.build(wp.cut(0, n / remaining_frac));
		// pkd.flatten(wp);
		// points_type wp2 = WP;
		// assert(parlay::sort(wp) == parlay::sort(wp2));
		// std::cout << "same points" << "\n";
	}

	return;
}

template <typename Point, typename Tree, bool test_time = true,
	  bool Serial = false>
void batch_insert(Tree &pkd, parlay::sequence<Point> const &WP,
		  parlay::sequence<Point> const &WI, int const &rounds,
		  double ratio = 1.0)
{
	using points_type = typename Tree::points_type;
	using box_type = typename Tree::box_type;
	points_type wp = points_type::uninitialized(WP.size());
	points_type wi = points_type::uninitialized(WI.size());
	double loop_late = rounds > 1 ? 0.01 : -100;

	// NOTE: build the tree by type
	auto build_tree_by_type = [&]() {
		if constexpr (psi::is_kd_tree<Tree> || psi::is_p_tree<Tree>) {
			parlay::copy(WP, wp), parlay::copy(WI, wi);
			pkd.build(parlay::make_slice(wp));
		} else if constexpr (psi::is_orth_tree<Tree>) {
			parlay::copy(WP, wp), parlay::copy(WI, wi);
			auto box1 = Tree::get_box(parlay::make_slice(wp));
			auto box2 = Tree::get_box(wi.cut(
				0, static_cast<size_t>(wi.size() * ratio)));
			box_type box = Tree::get_box(box1, box2);
			// std::cout << box1.first << ' ' << box1.second;
			// std::cout << box2.first << ' ' << box2.second;
			// std::cout << box.first << ' ' << box.second <<
			// std::endl;
			pkd.build(parlay::make_slice(wp), box);
		} else {
			parlay::copy(WP, wp), parlay::copy(WI, wi);
			pkd.build(parlay::make_slice(wp));
		}
	};

	if constexpr (test_time) { // NOTE: clean and measure time
		pkd.delete_tree();
		double aveInsert = time_loop(
			rounds, loop_late, [&]() { build_tree_by_type(); },
			[&]() {
				pkd.batch_insert(
					wi.cut(0, static_cast<size_t>(
							  wi.size() * ratio)));
			},
			[&]() { pkd.delete_tree(); });
		std::cout << aveInsert << " " << std::flush;
		build_tree<Point, Tree, false>(WP, rounds, pkd);
	} else { // NOTE: insert the points from previous tree
		pkd.delete_tree();
		build_tree_by_type();
		parlay::copy(WI, wi);
		pkd.batch_insert(
			wi.cut(0, static_cast<size_t>(wi.size() * ratio)));
		std::cout << "finish" << std::endl;
	}

	return;
}

template <typename Point, typename Tree, bool test_time = true,
	  bool Serial = false>
void batch_delete(Tree &pkd, parlay::sequence<Point> const &WP,
		  parlay::sequence<Point> const &WI, int const &rounds,
		  double ratio = 1.0)
{
	using points_type = typename Tree::points_type;
	points_type wp = points_type::uninitialized(WP.size());
	points_type wi = points_type::uninitialized(WP.size());
	size_t batchSize = static_cast<size_t>(WP.size() * ratio);
	double loop_late = rounds > 1 ? 0.01 : -100;

	if constexpr (test_time) {
		pkd.delete_tree();
		double aveDelete = time_loop(
			rounds, loop_late,
			[&]() {
				build_tree<Point, Tree, false>(WP, rounds, pkd);
				parlay::copy(WI, wi);
			},
			[&]() { pkd.batch_delete(wi.cut(0, batchSize)); },
			[&]() { pkd.delete_tree(); });
		std::cout << aveDelete << " " << std::flush;
		build_tree<Point, Tree, false>(WP, rounds, pkd);
	} else {
		parlay::copy(WI, wi);
		pkd.batch_delete(wi.cut(0, batchSize));
	}

	return;
}

template <typename Point, typename Tree, bool test_time = true>
void batch_diff(Tree &pkd, parlay::sequence<Point> const &WP, int const &rounds,
		double const total_ratio = 1.0,
		double const overlap_ratio = 0.5)
{
	using points_type = typename Tree::points_type;
	size_t total_batch_size = static_cast<size_t>(WP.size() * total_ratio);
	size_t overlap_size =
		static_cast<size_t>(total_batch_size * overlap_ratio);

	points_type WI =
		parlay::tabulate(total_batch_size, [&](size_t i) -> Point {
			if (i < overlap_size) {
				return WP[i];
			} else {
				Point p = WP[i];
				if constexpr (is_aug_point<Point>) {
					p.aug.id = -1 * p.aug.id - 1;
				} else {
					std::transform(
						p.pnt.begin(), p.pnt.end(),
						p.pnt.begin(),
						std::negate<coord_type>());
				}
				return p;
			}
		});

	points_type wp = points_type::uninitialized(WP.size());
	points_type wi = points_type::uninitialized(WI.size());

	if constexpr (test_time) {
		pkd.delete_tree();
		double aveDelete = time_loop(
			rounds, 0.01,
			[&]() {
				build_tree<Point, Tree, false>(WP, rounds, pkd);
				parlay::copy(WI, wi);
			},
			[&]() { pkd.batch_diff(wi.cut(0, total_batch_size)); },
			[&]() { pkd.delete_tree(); });
		std::cout << aveDelete << " " << std::flush;
		build_tree<Point, Tree, false>(WP, rounds, pkd);
	} else {
		parlay::copy(WI, wi);
		pkd.batch_diff(wi.cut(0, total_batch_size));
	}

	return;
}

struct step_update_logger {
	int id;
	double t;

	friend std::ostream &operator<<(std::ostream &os,
					step_update_logger const &log)
	{
		os << "(" << log.id << ", " << log.t << ")";
		return os;
	}
};

template <typename Point, typename Tree, bool Insert>
void batch_insert_by_step(Tree &pkd, parlay::sequence<Point> const &WP,
			  int const rounds, double const insert_ratio,
			  int const remain_divide_ratio = 2)
{
	using points_type = typename Tree::points_type;
	using box_type = typename Tree::box_type;
	points_type wp = points_type::uninitialized(WP.size());
	size_t n = wp.size();
	size_t step = static_cast<size_t>(insert_ratio * n);
	size_t slice_num = n / step;
	parlay::sequence<parlay::sequence<double>> time_table(
		rounds + 2, parlay::sequence<double>(slice_num, 0.0));
	parlay::sequence<step_update_logger> log_time(slice_num);
	int id_cnt = 0;
	for (auto &i : log_time) {
		i = {id_cnt++, 0.0};
	}
	size_t round_cnt = 0;

	pkd.delete_tree();

	// NOTE: build the tree by type
	auto prepare_build = [&]() {
		pkd.delete_tree();
		if constexpr (psi::is_kd_tree<Tree> || psi::is_p_tree<Tree>) {
			parlay::copy(WP, wp);
		} else if constexpr (psi::is_orth_tree<Tree>) {
			parlay::copy(WP, wp);
			auto box = Tree::get_box(wp.cut(0, n));
			pkd.set_bounding_box(box);
		} else {
			parlay::copy(WP, wp);
		}
	};

	auto incre_build = [&](int rounds_num_cutoff) {
		parlay::internal::timer t;
		size_t l = 0, r = 0;
		size_t cnt = 0;
		while (l < n) {
			if (cnt >= rounds_num_cutoff) {
				break;
			}

			r = std::min(l + step, n);
			// std::cout << l << " " << r << " " << n << std::endl;
			pkd.batch_insert(parlay::make_slice(wp.begin() + l,
							    wp.begin() + r));
			l = r;
			time_table[round_cnt][cnt++] += t.next_time();
		}
		round_cnt++;
	};

	double loop_late = rounds > 1 ? 0.01 : -100;
	double ave_time = time_loop(
		rounds, loop_late, [&]() { prepare_build(); },
		[&]() { incre_build(slice_num + 1); },
		[&]() { pkd.delete_tree(); });

	// begin count the time
	// if (rounds != 1 && round_cnt - 1 != rounds) {
	//   throw std::runtime_error("rounds not match!");
	// }
	for (int i = 1; i <= rounds; i++) {
		for (int j = 0; j < slice_num; j++) {
			log_time[j].t += time_table[i][j];
		}
	}
	for (int j = 0; j < slice_num; j++) {
		log_time[j].t /= rounds;
	}

	puts("# Insert");
	std::cout << "## " << insert_ratio << std::endl;
	std::sort(log_time.begin(), log_time.end(),
		  [](auto const &a, auto const &b) { return a.t < b.t; });

	// Calculate average time from log_time
	double total_time = 0.0;
	for (auto const &log : log_time) {
		total_time += log.t;
	}
	double average_time = total_time / log_time.size();

	std::cout << "median: " << log_time[slice_num / 2]
		  << "-> min: " << *log_time.begin()
		  << "-> max: " << *log_time.rbegin() << "-> tot: " << ave_time
		  << "-> avg: " << average_time << std::endl;

	// WARN: restore status
	prepare_build();
	incre_build(slice_num / remain_divide_ratio);
	// auto original_box = Tree::get_box(wp.cut(0, n /
	// remain_divide_ratio)); auto tree_box = pkd.get_root_box(); std::cout
	// << original_box.first << ' ' << original_box.second << std::endl;
	// std::cout << tree_box.first << ' ' << tree_box.second << std::endl;
	// auto root = pkd.cpam_aug_map_; std::cout << "root size: "
	// << root.size() << std::endl;

	return;
}

template <typename Point, typename Tree, bool Insert>
void batch_delete_by_step(Tree &pkd, parlay::sequence<Point> const &WP,
			  int const rounds, double const insert_ratio,
			  size_t const remain_divide_ratio = 2)
{
	using points_type = typename Tree::points_type;
	using box_type = typename Tree::box_type;
	points_type wp = points_type::uninitialized(WP.size());
	size_t n = wp.size();
	size_t step = static_cast<size_t>(insert_ratio * n);
	size_t slice_num = n / step;
	// std::cout << "n: " << n << " step: " << step << " slice_num: " <<
	// slice_num
	// << std::endl;
	parlay::sequence<parlay::sequence<double>> time_table(
		rounds + 2, parlay::sequence<double>(slice_num, 0.0));
	parlay::sequence<step_update_logger> log_time(slice_num);
	int id_cnt = 0;
	for (auto &i : log_time) {
		i = {id_cnt++, 0.0};
	}
	size_t round_cnt = 0;

	pkd.delete_tree();

	// NOTE: build the tree by type
	auto build_tree_by_type = [&]() {
		parlay::copy(WP, wp);

		if constexpr (psi::is_kd_tree<Tree> || psi::is_p_tree<Tree>) {
			pkd.build(parlay::make_slice(wp));
		} else if constexpr (psi::is_orth_tree<Tree>) {
			auto box = Tree::get_box(wp.cut(0, n));
			pkd.build(parlay::make_slice(wp), box);
		} else {
			pkd.build(parlay::make_slice(wp));
		}

		parlay::copy(WP, wp);
	};

	auto incre_delete = [&](int rounds_num_cutoff, auto wp_slice) {
		parlay::internal::timer t;
		size_t l = 0, r = 0;
		size_t cnt = 0;
		while (l < n) {
			if (cnt >= rounds_num_cutoff) {
				break;
			}

			r = std::min(l + step, n);
			// WARN: r may exceeds the right bounds with offset
			pkd.batch_delete(parlay::make_slice(
				wp_slice.begin() + l, wp_slice.begin() + r));
			l = r;
			time_table[round_cnt][cnt++] += t.next_time();
		}
		round_cnt++;
	};

	double loop_late = rounds > 1 ? 0.01 : -100;
	double ave_time = time_loop(
		rounds, loop_late, [&]() { build_tree_by_type(); },
		[&]() { incre_delete(slice_num, wp.cut(0, slice_num * step)); },
		[&]() { pkd.delete_tree(); });

	// if (rounds != 1 && round_cnt - 1 != rounds) {
	//   throw std::runtime_error("rounds not match!");
	// }
	for (int i = 1; i <= rounds; i++) {
		for (int j = 0; j < slice_num; j++) {
			log_time[j].t += time_table[i][j];
		}
	}
	for (int j = 0; j < slice_num; j++) {
		log_time[j].t /= rounds;
	}

	puts("# Delete");
	std::cout << "## " << insert_ratio << std::endl;
	// puts("");
	// std::cout << ave_time << " " << std::flush;
	std::sort(log_time.begin(), log_time.end(),
		  [](auto const &a, auto const &b) { return a.t < b.t; });
	// Calculate average time from log_time
	double total_time = 0.0;
	for (auto const &log : log_time) {
		total_time += log.t;
	}
	double average_time = total_time / log_time.size();

	std::cout << "median: " << log_time[slice_num / 2]
		  << "-> min: " << *log_time.begin()
		  << "-> max: " << *log_time.rbegin() << "-> tot: " << ave_time
		  << "-> avg: " << average_time << std::endl;

	// WARN: restore status
	pkd.delete_tree();
	build_tree_by_type();
	incre_delete(slice_num / remain_divide_ratio,
		     wp.cut(slice_num / remain_divide_ratio * step, n));

	return;
}

template <typename Point, typename Tree, bool printHeight = 0,
	  bool printVisNode = 1>
void queryKNN([[maybe_unused]] uint_fast8_t const &Dim,
	      parlay::sequence<Point> const &WP, int const &rounds, Tree &pkd,
	      scalar_type *kdknn, int const K, bool const flattenTreeTag)
{
	using points_type = typename Tree::points_type;
	using coord_type = typename Point::coord_type;
	using dis_type = typename Point::dis_type;
	using nn_pair = std::pair<std::reference_wrapper<Point>, dis_type>;
	using leaf_type = typename Tree::leaf_type;
	using interior_type = typename Tree::interior_type;
	// using nn_pair = std::pair<Point, dis_type>;
	size_t n = WP.size();
	// int LEAVE_WRAP = 32;
	double loopLate = rounds > 1 ? 0.01 : -0.1;
	auto *kd_parallel_root = pkd.get_root();

	points_type wp = WP;
	// points_type wp = parlay::random_shuffle(WP);

	parlay::sequence<nn_pair> out(
		K * n, nn_pair(std::ref(wp[0]), static_cast<dis_type>(0)));
	// parlay::sequence<nn_pair> out(K * n);
	parlay::sequence<bounded_queue<Point, nn_pair>> bq =
		parlay::sequence<bounded_queue<Point, nn_pair>>::uninitialized(
			n);
	parlay::parallel_for(0, n, [&](size_t i) {
		bq[i].resize(out.cut(i * K, i * K + K));
	});
	parlay::sequence<size_t> vis_leaf(n), vis_inter(n), gen_box(n),
		check_box(n), skip_box(n);

	double aveQuery = time_loop(
		rounds, loopLate,
		[&]() {
			parlay::parallel_for(0, n,
					     [&](size_t i) { bq[i].reset(); });
		},
		[&]() {
			// if (!flattenTreeTag) {  // WARN: Need ensure
			// pkd.size() == wp.size()
			//   pkd.flatten(parlay::make_slice(wp));
			// }
			parlay::parallel_for(0, n, [&](size_t i) {
				// for (size_t i = 0; i < n; i++) {
				/* The baseline trees still take the root; the
				 * psi trees dropped it, every caller was
				 * passing get_root(). */
				auto knn_log = [&]() {
					if constexpr (psi::is_kd_tree<Tree> ||
						      psi::is_orth_tree<Tree> ||
						      psi::is_p_tree<Tree>) {
						return pkd.knn(wp[i], bq[i]);
					} else {
						return pkd.knn(kd_parallel_root,
							       wp[i], bq[i]);
					}
				}();
				auto [vis_leaf_num, vis_inter_num, gen_box_num,
				      check_box_num, skip_box_num] = knn_log;
				kdknn[i] = bq[i].top().second;
				vis_leaf[i] = vis_leaf_num;
				vis_inter[i] = vis_inter_num;
				gen_box[i] = gen_box_num;
				check_box[i] = check_box_num;
				skip_box[i] = skip_box_num;
				// }
			});
		},
		[&]() {});

	std::cout << aveQuery << " " << std::flush;
	if (printHeight) {
		// WARN: change when using multi-node
		size_t max_deep = 0;
		std::cout << pkd.template get_max_tree_depth<leaf_type,
							     interior_type>(
				     pkd.get_root(), max_deep)
			  << " "
			  << pkd.template get_ave_tree_height<leaf_type,
							      interior_type>()
			  << " " << std::flush;
	}
	if (printVisNode) {
		std::cout << static_cast<double>(
				     parlay::reduce(vis_leaf.cut(0, n))) /
				     static_cast<double>(n)
			  << " " << std::flush;
		std::cout << static_cast<double>(
				     parlay::reduce(vis_inter.cut(0, n))) /
				     static_cast<double>(n)
			  << " " << std::flush;
		std::cout << static_cast<double>(
				     parlay::reduce(gen_box.cut(0, n))) /
				     static_cast<double>(n)
			  << " " << std::flush;
		std::cout << static_cast<double>(
				     parlay::reduce(check_box.cut(0, n))) /
				     static_cast<double>(n)
			  << " " << std::flush;
		std::cout << static_cast<double>(
				     parlay::reduce(skip_box.cut(0, n))) /
				     static_cast<double>(n)
			  << " " << std::flush;
	}

	return;
}

// NOTE: run range count and check the correct
template <typename Point, typename Tree>
void range_count(parlay::sequence<Point> const &wp, Tree &pkd,
		 int const &rounds, int rec_num, int rec_type, int const DIM)
{
	// using points_type = typename Tree::points_type;
	// using box_type = typename Tree::box_type;

	auto [query_box_seq, max_size] =
		gen_rectangles<Point, Tree, false>(rec_num, rec_type, wp, DIM);
	parlay::sequence<size_t> kdknn(rec_num, 0);

	double aveCount = time_loop(
		rounds, 0.01, [&]() {},
		[&]() {
			parlay::parallel_for(0, rec_num, [&](size_t i) {
				auto [size, logger] =
					pkd.range_count(query_box_seq[i].first);
				kdknn[i] = size;
			});
		},
		[&]() {});

	// NOTE: verify the solutions
	std::cout << "check range count: " << rec_num << " " << rec_type
		  << std::endl;
	for (int i = 0; i < rec_num; i++) {
		if (kdknn[i] != query_box_seq[i].second) {
			std::cout << kdknn[i] << " " << query_box_seq[i].second
				  << " " << query_box_seq[i].first.first << " "
				  << query_box_seq[i].first.second << std::endl;
		}
		assert(std::cmp_equal(kdknn[i], query_box_seq[i].second));
	}

	std::cout << aveCount << "\n" << std::flush;

	return;
}

template <typename Point, typename Tree>
void rangeCountRadius(parlay::sequence<Point> const &wp, Tree &pkd,
		      scalar_type *kdknn, int const &rounds,
		      int const &queryNum)
{
	// using points_type = typename Tree::points_type;
	// using node = typename Tree::node;
	using box_type = typename Tree::box_type;
	using circle = typename Tree::circle;

	int n = wp.size();

	double aveCount = time_loop(
		rounds, 0.01, [&]() {},
		[&]() {
			parlay::parallel_for(0, queryNum, [&](size_t i) {
				box_type query_box_seq = pkd.get_box(
					box_type(wp[i], wp[i]),
					box_type(wp[(i + n / 2) % n],
						 wp[(i + n / 2) % n]));
				auto d = Tree::p2p_distance(wp[i],
							    wp[(i + n / 2) % n],
							    wp[i].get_dim());
				d = static_cast<coord_type>(std::sqrt(d));
				circle cl = circle(wp[i], d);
				kdknn[i] = pkd.range_count(cl);
			});
		},
		[&]() {});

	std::cout << aveCount << " " << std::flush;

	return;
}

// NOTE: run range query and check the correct
template <typename Point, typename Tree>
void range_query(parlay::sequence<Point> const &wp, Tree &pkd,
		 int const &rounds, int const rec_num, int const rec_type,
		 int const DIM)
{
	using points_type = typename Tree::points_type;

	auto [query_box_seq, max_size] =
		gen_rectangles<Point, Tree, true>(rec_num, rec_type, wp, DIM);
	auto [offset, tot_size] = parlay::scan(
		parlay::delayed_tabulate(
			rec_num,
			[&](size_t i) -> size_t {
				return query_box_seq[i].second.size();
			}),
		parlay::addm<size_t>());
	offset.push_back(tot_size);
	points_type out(tot_size);
	parlay::sequence<size_t> kdknn(rec_num, 0);
	parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
		gen_box(rec_num), full_box(rec_num), skip_box(rec_num);

	double aveQuery = time_loop(
		rounds, 0.01, [&]() {},
		[&]() {
			parlay::parallel_for(0, rec_num, [&](size_t i) {
				auto [size, logger] = pkd.range_query(
					query_box_seq[i].first,
					out.cut(offset[i], offset[i + 1]));
				kdknn[i] = size;
				vis_leaf[i] = logger.vis_leaf_num;
				vis_inter[i] = logger.vis_interior_num;
				gen_box[i] = logger.generate_box_num;
				full_box[i] = logger.full_box_num;
				skip_box[i] = logger.skip_box_num;
			});
		},
		[&]() {});

	std::cout << "check range query: " << rec_num << " " << rec_type << " "
		  << max_size << std::endl;
	for (int i = 0; i < rec_num; i++) {
		if (kdknn[i] != query_box_seq[i].second.size()) {
			std::cout << kdknn[i] << " "
				  << query_box_seq[i].second.size() << " "
				  << query_box_seq[i].first.first
				  << query_box_seq[i].first.second << std::endl;
			std::cout << vis_leaf[i] << " " << gen_box[i] << " "
				  << full_box[i] << " " << skip_box[i]
				  << std::endl;
		}
		assert(std::cmp_equal(kdknn[i],
				      query_box_seq[i].second.size()));
		// std::cout << kdknn[i] << " " <<
		// query_box_seq[i].second.size() << " "
		//     << query_box_seq[i].first.first <<
		//     query_box_seq[i].first.second
		//     << std::endl;
		parlay::sort_inplace(out.cut(offset[i], offset[i + 1]),
				     [&](auto const &a, auto const &b) {
					     return a.aug.id < b.aug.id;
				     });
		parlay::sort_inplace(query_box_seq[i].second,
				     [&](auto const &a, auto const &b) {
					     return a.aug.id < b.aug.id;
				     });
		for (size_t j = 0; j < query_box_seq[i].second.size(); j++) {
			if (out[offset[i] + j] !=
			    query_box_seq[i].second.at(j)) {
				std::cout << "wrong"
					  << query_box_seq[i].first.first
					  << query_box_seq[i].first.second
					  << std::endl;
				std::cout << out[offset[i] + j] << " "
					  << query_box_seq[i].second.at(j)
					  << std::endl;
			}

			if constexpr (is_kd_tree<Tree> ||
				      is_orth_tree<Tree>) { // TODO: fix this by
							    // enable kdtree
				// handling duplicates
				// by id
				// assert(out[offset[i] +
				// j].same_dimension(query_box_seq[i].second.at(j)));
				assert(out[offset[i] + j] ==
				       query_box_seq[i].second.at(j));
			} else if constexpr (is_p_tree<Tree>) {
				assert(out[offset[i] + j] ==
				       query_box_seq[i].second.at(j));
			}

			// if (out[i * step + j] !=
			// query_box_seq[i].second.at(j)) std::cout << "wrong
			// "; std::cout << out[j] << " " <<
			// query_box_seq[i].second.at(j) << std::endl;
		}
	}

	std::cout << aveQuery << "\n" << std::flush;
	return;
}

//* test range count for fix rectangle
template <typename Point, typename Tree>
void rangeCountFix(Tree &pkd, scalar_type *kdknn, int const &rounds,
		   int rec_type, int rec_num, int DIM,
		   auto const &query_box_seq, auto max_size)
{
	// using Tree = Tree;
	// using points_type = typename Tree::points_type;
	// using box_type = typename Tree::box_type;

	// int n = WP.size();

	// auto [query_box_seq, max_size] =
	//     gen_rectangles<Point, Tree, false>(rec_num, rec_type, WP, DIM);
	parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
		gen_box(rec_num), full_box(rec_num), skip_box(rec_num);

	double aveCount = time_loop(
		rounds, 0.01, [&]() {},
		[&]() {
			parlay::parallel_for(0, rec_num, [&](size_t i) {
				auto [size, logger] =
					pkd.range_count(query_box_seq[i].first);

				kdknn[i] = size;
				vis_leaf[i] = logger.vis_leaf_num;
				vis_inter[i] = logger.vis_interior_num;
				gen_box[i] = logger.generate_box_num;
				full_box[i] = logger.full_box_num;
				skip_box[i] = logger.skip_box_num;
			});
			// for (int i = 0; i < rec_num; i++) {
			//     kdknn[i] =
			//     pkd.range_count(query_box_seq[i].first);
			// }
		},
		[&]() {});

	std::cout << aveCount << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(vis_leaf.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(vis_inter.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(gen_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(full_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(skip_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;

	return;
}
//
// template<typename Point>
// void rangeCountFixWithLog(const parlay::sequence<Point>& WP,
// base_tree<Point>& pkd, scalar_type* kdknn, const int& rounds,
//                           int rec_type, int rec_num, int DIM) {
//     using Tree = base_tree<Point>;
//     using points_type = typename Tree::points_type;
//     using node = typename Tree::node;
//     using box_type = typename Tree::box_type;
//
//     int n = WP.size();
//
//     auto [query_box_seq, max_size] = gen_rectangles(rec_num, rec_type, WP,
//     DIM); parlay::sequence<size_t> visLeafNum(rec_num, 0),
//     visInterNum(rec_num, 0); parlay::internal::timer t; for (int i = 0; i <
//     rec_num; i++) {
//         double aveQuery = time_loop(
//             rounds, -1.0,
//             [&]() {
//                 visInterNum[i] = 0;
//                 visLeafNum[i] = 0;
//             },
//             [&]() { kdknn[i] = pkd.range_count(query_box_seq[i].first,
//             visLeafNum[i], visInterNum[i]); }, [&]() {});
//         if (query_box_seq[i].second != kdknn[i]) std::cout << "wrong" <<
//         std::endl; std::cout << query_box_seq[i].second << " " <<
//         std::scientific << aveQuery
//         << std::endl;
//     }
//
//     return;
// }
//
//* test range query for fix rectangle
template <typename Point, typename Tree>
void rangeQueryFix(Tree &pkd, scalar_type *kdknn, int const &rounds,
		   parlay::sequence<Point> &out, int rec_type, int rec_num,
		   int DIM, auto const &query_box_seq, auto max_size)
{
	// auto [query_box_seq, max_size] =
	//     gen_rectangles<Point, Tree, false>(rec_num, rec_type, WP, DIM);
	parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
		gen_box(rec_num), full_box(rec_num), skip_box(rec_num);
	auto [offset, tot_size] = parlay::scan(
		parlay::delayed_tabulate(rec_num,
					 [&](size_t i) -> size_t {
						 return query_box_seq[i].second;
					 }),
		parlay::addm<size_t>());
	offset.push_back(tot_size);
	// std::cout << "range query: " << rec_num << " " << rec_type << " " <<
	// tot_size
	//           << std::endl;
	out.resize(tot_size);

	// int n = WP.size();
	// size_t step = out.size() / rec_num;
	// using ref_t = std::reference_wrapper<Point>;
	// parlay::sequence<ref_t> out_ref( out.size(), std::ref( out[0] ) );

	double aveQuery = time_loop(
		rounds, 0.01, [&]() {},
		[&]() {
			parlay::parallel_for(0, rec_num, [&](size_t i) {
				auto [size, logger] = pkd.range_query(
					query_box_seq[i].first,
					out.cut(offset[i], offset[i + 1]));

				kdknn[i] = size;
				vis_leaf[i] = logger.vis_leaf_num;
				vis_inter[i] = logger.vis_interior_num;
				gen_box[i] = logger.generate_box_num;
				full_box[i] = logger.full_box_num;
				skip_box[i] = logger.skip_box_num;
			});
		},
		[&]() {});

	std::cout << aveQuery << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(vis_leaf.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(vis_inter.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(gen_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(full_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	std::cout << static_cast<double>(
			     parlay::reduce(skip_box.cut(0, rec_num))) /
			     static_cast<double>(rec_num)
		  << " " << std::flush;
	return;
}

template <typename Point, typename Tree>
void range_query_serial_with_log(Tree &pkd, scalar_type *kdknn,
				 int const &rounds,
				 parlay::sequence<Point> &out, int rec_type,
				 int rec_num, int DIM,
				 auto const &query_box_seq, auto max_size)
{
	auto [offset, tot_size] = parlay::scan(
		parlay::delayed_tabulate(rec_num,
					 [&](size_t i) -> size_t {
						 return query_box_seq[i].second;
					 }),
		parlay::addm<size_t>());
	offset.push_back(tot_size);
	// using ref_t = std::reference_wrapper<point>;
	// parlay::sequence<ref_t> out_ref( out.size(), std::ref( out[0] ) );
	out.resize(tot_size);

	for (int i = 0; i < rec_num; i++) {
		parlay::internal::timer t;
		t.reset(), t.start();
		auto [size, logger] =
			pkd.range_query(query_box_seq[i].first,
					out.cut(offset[i], offset[i + 1]));
		t.stop();
		// double ave_query = time_loop(
		//     rounds, -1.0, [&]() {},
		//     [&]() {
		//       auto [size, logger] = pkd.range_query(
		//           query_box_seq[i].first, out.cut(offset[i], offset[i
		//           + 1]));
		//     },
		//     [&]() {});
		std::cout << rec_type << " " << query_box_seq[i].second << " "
			  << std::scientific << t.total_time() << std::endl;
	}

	return;
}

template <typename T>
class counter_iterator
{
private:
	struct accept_any {
		template <typename U>
		accept_any &operator=(U const &)
		{
			return *this;
		}
	};

public:
	typedef std::output_iterator_tag iterator_category;

	counter_iterator(T &counter) : counter(counter)
	{
	}
	counter_iterator &operator=(counter_iterator const &other)
	{
		if (this != &other) { // Check for self-assignment
			counter = other.counter;
		}
		return *this;
	}
	counter_iterator(counter_iterator const &other) : counter(other.counter)
	{
	}

	bool operator==(counter_iterator const &rhs) const
	{
		return counter == rhs.counter;
	}
	bool operator!=(counter_iterator const &rhs) const
	{
		return counter != rhs.counter;
	}

	accept_any operator*() const
	{
		++counter.get();
		return {};
	}

	counter_iterator &operator++()
	{ // ++a
		return *this;
	}
	counter_iterator operator++(int)
	{ // a++
		return *this;
	}

protected:
	std::reference_wrapper<T> counter;
};

//*---------- generate points_type within a 0-box_size --------------------
template <typename Point>
void generate_random_points(parlay::sequence<Point> &wp, coord_type _box_size,
			    long n, int Dim)
{
	coord_type box_size = _box_size;

	std::mt19937 gen_mt(generator_seed());
	std::uniform_int_distribution<int> distrib(1, box_size);

	parlay::random_generator gen(distrib(gen_mt));
	std::uniform_int_distribution<int> dis(0, box_size);

	wp.resize(n);
	// generate n random points_type in a cube
	parlay::parallel_for(
		0, n,
		[&](long i) {
			auto r = gen[i];
			for (int j = 0; j < Dim; j++) {
				wp[i][j] = dis(r);
			}
		},
		1000);
	return;
}

template <typename Point>
std::pair<size_t, int> read_points(char const *iFile,
				   parlay::sequence<Point> &wp,
				   [[maybe_unused]] int id_offset = 0)
{
	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	static coords_type a_sample_point;
	parlay::sequence<char> S = readStringFromFile(iFile);
	parlay::sequence<char *> W = stringToWords(S);
	size_t N = std::stoul(W[0], nullptr, 10);
	int Dim = atoi(W[1]);
	assert(N > 0 && Dim >= 1);

	auto pts = W.cut(2, W.size());
	assert(pts.size() % Dim == 0);
	size_t n = pts.size() / Dim;
	auto a = parlay::tabulate(Dim * n, [&](size_t i) -> coord_type {
		if constexpr (std::is_integral_v<coord_type>)
			return std::stol(pts[i]);
		else if (std::is_floating_point_v<coord_type>)
			return std::stod(pts[i]);
	});
	wp.resize(N);
	parlay::parallel_for(0, n, [&](size_t i) {
		for (int j = 0; j < Dim; j++) {
			wp[i][j] = a[i * Dim + j];
			if constexpr (std::is_same_v<
					      Point,
					      basic_point<
						      coord_type,
						      a_sample_point.size()>>) {
				;
			} else if constexpr (
				std::is_same_v<Point,
					       typename ZD::geobase::Point> ||
				std::is_same_v<Point,
					       typename ZD3D::geobase::Point>) {
				wp[i].id = i + id_offset;
			} else {
				wp[i].aug.id = i + id_offset;
			}
		}
	});
	return std::make_pair(N, Dim);
}

template <typename TreeWrapper>
void print_tree_param()
{
	std::cout << "Tree: " << TreeWrapper::tree_type::get_tree_name() << "; "
		  << "AugType: " << TreeWrapper::tree_type::check_has_box()
		  << "; "
		  << "Split: " << TreeWrapper::SplitRule::get_split_name()
		  << "; "
		  << "BDO: " << TreeWrapper::tree_type::get_build_depth_once()
		  << "; "
		  << "Inba: " << TreeWrapper::tree_type::get_imbalance_ratio()
		  << "; ";

	if constexpr (std::is_integral_v<
			      typename TreeWrapper::Point::coord_type>) {
		std::cout << "Coord: integer"
			  << "; ";
	} else if (std::is_floating_point_v<
			   typename TreeWrapper::Point::coord_type>) {
		std::cout << "Coord: float"
			  << "; ";
	}
	std::cout << "\n" << std::flush;
	return;
}

// NOTE: default test functions for all custom tree
static auto constexpr default_test_func = []<class TreeDesc, typename Point>(
						  int const &num_dims,
						  parlay::sequence<Point> const
							  &wp,
						  parlay::sequence<Point> const
							  &wi,
						  size_t const &N, int const &K,
						  int const &kRounds,
						  string const &kInsertFile,
						  int const &kTag,
						  int const &kQueryType,
						  int const kSummary) {
	using Tree = TreeDesc::tree_type;
	using points_type = typename Tree::points_type;

	Tree tree;
	constexpr bool test_time = true;

	// std::cout << "Called build" << std::endl;
	// build_tree<Point, Tree, test_time, 2>(wp, kRounds, tree);
	// std::cout << "build Finished" << std::endl;

	// NOTE: batch insert
	if (kTag & (1 << 0)) {
		if (kSummary) {
			parlay::sequence<double> const ratios = {0.0001, 0.001,
								 0.01, 0.1};
			for (size_t i = 0; i < ratios.size(); i++) {
				batch_insert<Point, Tree, test_time>(
					tree, wp, wi, kRounds, ratios[i]);
			}
		} else {
			batch_insert<Point, Tree, test_time>(tree, wp, wi,
							     kRounds, 0.1);
		}
	}

	// NOTE: batch delete
	if (kTag & (1 << 1)) {
		if (kSummary) {
			parlay::sequence<double> const ratios = {0.0001, 0.001,
								 0.01, 0.1};
			for (size_t i = 0; i < ratios.size(); i++) {
				batch_delete<Point, Tree, test_time>(
					tree, wp, wp, kRounds, ratios[i]);
			}
		} else {
			batch_delete<Point, Tree, test_time>(
				tree, wp, wp, kRounds, batch_insert_ratio);
		}
	}

	scalar_type *kdknn = nullptr;
	auto run_batch_knn = [&](points_type const &query_pts, int kth) {
		kdknn = new scalar_type[query_pts.size()];
		queryKNN<Point>(num_dims, query_pts, kRounds, tree, kdknn, kth,
				true);
		delete[] kdknn;
	};

	auto generate_query_box = [&](int rec_num, int rec_total_type,
				      points_type const &wp) {
		// NOTE: generate rectangles for the first half of the points
		parlay::sequence<parlay::sequence<
			std::pair<typename Tree::box_type, size_t>>>
			query_box_seq(rec_total_type);
		parlay::sequence<size_t> query_max_size(rec_total_type);
		for (int i = 0; i < rec_total_type; i++) {
			auto [query_box, max_size] =
				gen_rectangles<Point, Tree, false, true>(
					rec_num, i, wp, num_dims);
			query_box_seq[i] = query_box;
			query_max_size[i] = max_size;
		}
		return std::make_pair(query_box_seq, query_max_size);
	};

	auto incre_update_test_bundle = [&](auto const &query_box_seq,
					    auto const &query_max_size) {
		// NOTE: knn query
		{
			int k[3] = {1, 10, 100};

			std::cout << "in-dis-skewed knn time: ";
			size_t batch_size = static_cast<size_t>(
				wp.size() * batch_query_ratio);
			for (int i = 0; i < 3; i++) {
				run_batch_knn(wp.subseq(0, batch_size), k[i]);
			}
			puts("");

			std ::cout << "out-dis-skewed knn time: ";
			for (int i = 0; i < 3; i++) {
				run_batch_knn(wp.subseq(wp.size() - batch_size,
							wp.size()),
					      k[i]);
			}
			puts("");

			// NOTE: sample points within the whole input datasets
			auto query_pts = parlay::pack(
				wp,
				parlay::tabulate(
					wp.size(), [&](size_t i) -> bool {
						return i % (wp.size() /
							    (batch_size * 2)) ==
						       0;
					}));

			std::cout << "in-dis-uniform knn time: ";
			for (int i = 0; i < 3; i++) {
				run_batch_knn(
					parlay::random_shuffle(query_pts.subseq(
						0, batch_size)),
					k[i]);
			}
			puts("");

			std::cout << "out-dis-uniform knn time: ";
			for (int i = 0; i < 3; i++) {
				run_batch_knn(
					parlay::random_shuffle(query_pts.subseq(
						batch_size, query_pts.size())),
					k[i]);
			}
			puts("");
		}

		// NOTE: range count
		{
			int rec_num = query_box_seq[0].size();
			kdknn = new scalar_type[rec_num];

			std::cout << "range count time: ";
			for (int i = 0; i < 3; i++) {
				rangeCountFix<Point>(tree, kdknn, kRounds, i,
						     rec_num, num_dims,
						     query_box_seq[i],
						     query_max_size[i]);
			}
			delete[] kdknn;
			puts("");
		}

		// NOTE: range query
		{
			int rec_num = query_box_seq[0].size();
			kdknn = new scalar_type[rec_num];

			std::cout << "range query time: ";
			for (int i = 0; i < 3; i++) {
				points_type out;
				rangeQueryFix<Point>(tree, kdknn, kRounds, out,
						     i, rec_num, num_dims,
						     query_box_seq[i],
						     query_max_size[i]);
			}
			delete[] kdknn;
			puts("");
		}
	};

	// NOTE: scalability
	if (kTag & (1 << 2)) {
		puts("");
		build_tree<Point, Tree, test_time, 0>(wp, kRounds, tree);
		batch_insert<Point, Tree, test_time>(tree, wp, wi, kRounds,
						     batch_insert_ratio);
		batch_delete<Point, Tree, test_time>(tree, wp, wp, kRounds,
						     batch_insert_ratio);
	}

	// NOTE: batch insert by step
	if (kTag & (1 << 3)) {
		puts("");
		build_tree<Point, Tree, test_time, 3>(wp, kRounds, tree, 2);

		auto [query_box_seq, query_max_size] = generate_query_box(
			range_query_num, 3, wp.subseq(0, wp.size() / 2));

		incre_update_test_bundle(query_box_seq, query_max_size);

		parlay::sequence<double> const ratios = {0.1, 0.01, 0.001,
							 0.0001};
		for (auto rat : ratios) {
			batch_insert_by_step<Point, Tree, true>(tree, wp,
								kRounds, rat);
			incre_update_test_bundle(query_box_seq, query_max_size);
		}
	}

	// NOTE: batch delete by step
	if (kTag & (1 << 4)) {
		puts("");
		auto [query_box_seq, query_max_size] = generate_query_box(
			range_query_num, 3, wp.subseq(0, wp.size() / 2));
		build_tree<Point, Tree, test_time, 3>(wp, kRounds, tree, 2);
		incre_update_test_bundle(query_box_seq, query_max_size);

		parlay::sequence<double> const ratios = {0.1, 0.01, 0.001,
							 0.0001};
		// parlay::sequence<double> const ratios = {0.001};
		for (auto rat : ratios) {
			batch_delete_by_step<Point, Tree, true>(tree, wp,
								kRounds, rat);
			incre_update_test_bundle(query_box_seq, query_max_size);
			// incre_update_test_bundle(wp.subseq(0, wp.size() /
			// 2));
		}
	}

	// real world
	if (kTag & (1 << 5)) {
		puts("");
		auto [query_box_seq, query_max_size] = generate_query_box(
			range_query_num, 3, wp.subseq(0, wp.size() / 2));

		build_tree<Point, Tree, test_time, 3>(wp, kRounds, tree, 2);
		incre_update_test_bundle(query_box_seq, query_max_size);

		batch_insert_by_step<Point, Tree, true>(tree, wp, kRounds,
							0.0001);
		incre_update_test_bundle(query_box_seq, query_max_size);

		batch_delete_by_step<Point, Tree, true>(tree, wp, kRounds,
							0.0001);
		incre_update_test_bundle(query_box_seq, query_max_size);
	}

	// range query with log
	if (kTag & (1 << 6)) {
		puts("");
		batch_insert_by_step<Point, Tree, true>(tree, wp, kRounds,
							0.0001);

		auto [query_box_seq, query_max_size] =
			generate_query_box(single_query_log_repeat_num, 3,
					   wp.subseq(0, wp.size() / 2));

		// NOTE: range query
		{
			int rec_num = single_query_log_repeat_num;
			kdknn = new scalar_type[rec_num];

			std::cout << "range query time: " << std::endl;
			for (int i = 0; i < 3; i++) {
				// std::cout << "range query time: " <<
				// std::endl;
				points_type out;
				// rangeQueryFix<Point>(tree, kdknn, kRounds,
				// out, i, rec_num, num_dims,
				//                      query_box_seq[i],
				//                      query_max_size[i]);
				range_query_serial_with_log<Point>(
					tree, kdknn, kRounds, out, i, rec_num,
					num_dims, query_box_seq[i],
					query_max_size[i]);
			}
			delete[] kdknn;
			puts("");
		}
	}

	// Batch updates
	if (kTag & (1 << 7)) {
		puts("");
		parlay::sequence<double> const ratios = {
			0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01,
			0.02,	0.05,	0.1,	0.2,   0.5,   1.0};
		// std::cout << "Insert: ";
		for (size_t i = 0; i < ratios.size(); i++) {
			batch_insert<Point, Tree, test_time>(
				tree, wp, wi, kRounds, ratios[i]);
		}
		// std::cout << std::endl;
		// std::cout << "Delete: ";
		for (size_t i = 0; i < ratios.size(); i++) {
			batch_delete<Point, Tree, test_time>(
				tree, wp, wp, kRounds, ratios[i]);
		}
		puts("");
	}
	// WARN: compress the kdnode to multi_node, should remove except for
	// exp if constexpr (is_kd_tree<Tree>) {
	//   tree.Compress2Multi();
	// }

	// batch_insert_by_step<Point, Tree, true>(tree, wp, kRounds,
	// 0.000000001);
	if (kQueryType & (1 << 0)) { // NOTE: knn
		size_t batch_size =
			static_cast<size_t>(wp.size() * batch_query_ratio);

		if (kSummary == 0) {
			int k[3] = {1, 10, 100};
			for (int i = 0; i < 3; i++) {
				run_batch_knn(wp.subseq(0, batch_size), k[i]);
			}
		} else { // test kSummary
			run_batch_knn(wp.subseq(0, batch_size), K);
		}
	}

	if (kQueryType & (1 << 1)) { // NOTE: range count
		auto [query_box_seq, query_max_size] =
			generate_query_box(range_query_num, 3, wp);

		if (!kSummary) {
			int recNum = range_query_num;
			kdknn = new scalar_type[recNum];

			// std::cout << std::endl;
			for (int i = 0; i < 3; i++) {
				rangeCountFix<Point>(tree, kdknn, kRounds, i,
						     recNum, num_dims,
						     query_box_seq[i],
						     query_max_size[i]);
			}

			delete[] kdknn;
		}
	}

	if (kQueryType & (1 << 2)) { // NOTE: range query
		if (kSummary == 0) {
			auto [query_box_seq, query_max_size] =
				generate_query_box(range_query_num, 3, wp);

			int recNum = range_query_num;
			kdknn = new scalar_type[recNum];

			for (int i = 0; i < 3; i++) {
				points_type out;
				rangeQueryFix<Point>(tree, kdknn, kRounds, out,
						     i, recNum, num_dims,
						     query_box_seq[i],
						     query_max_size[i]);
			}
			delete[] kdknn;
		} else if (kSummary == 1) { // NOTE: for kSummary
			auto [query_box_seq, query_max_size] =
				generate_query_box(summary_range_query_num, 3,
						   wp);

			kdknn = new scalar_type[summary_range_query_num];
			points_type out;
			rangeQueryFix<Point>(tree, kdknn, kRounds, out, 2,
					     summary_range_query_num, num_dims,
					     query_box_seq[2],
					     query_max_size[2]);
			delete[] kdknn;
		}
	}

	std::cout << "\n" << std::flush;

	tree.delete_tree();

	return;
};

class wrapper
{
public:
	// NOTE: Trees
	template <class PointType, class split_rule_type, class LeafAugType,
		  class InteriorAugType>
	struct kd_tree_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type = psi::kd_tree<psi::tree_traits<
			Point, SplitRule, LeafAugType, InteriorAugType>>;
	};

	template <class PointType, class split_rule_type, class LeafAugType,
		  class InteriorAugType>
	struct orth_tree_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type = psi::orth_tree<psi::orth_tree_traits<
			Point, SplitRule, LeafAugType, InteriorAugType>>;
	};

	template <class PointType, class split_rule_type>
	struct p_tree_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type =
			psi::p_tree<psi::tree_traits<Point, SplitRule>>;
	};

	template <class PointType, class split_rule_type>
	struct cpam_raw_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type =
			typename cpam_tree::cpam_raw<Point, SplitRule>;
	};

	/* zdtree wrapper */
	template <class PointType, class split_rule_type>
	struct zd_tree_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type = typename ZD::zdtree<Point, SplitRule>;
	};

	// For zdtree_3d
	template <class PointType, class split_rule_type>
	struct zd_tree_3d_wrapper {
		using Point = PointType;
		using SplitRule = split_rule_type;
		using tree_type = typename ZD3D::zdtree<Point, SplitRule>;
	};

	// NOTE: apply the dim and split rule
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
		friend std::ostream &operator<<(std::ostream &os,
						aug_id const &rhs)
		{
			os << rhs.id;
			return os;
		}
	};

	// For the spatial filling curve, we use the aug_id_code to
	// ensure the id is unique and the code is used to determine the
	// order of the points in the tree.
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
			// return code == rhs.code && id == rhs.id;
			// WARN: code is not important, we only need to ensure
			// the id
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

	// NOTE: driven functions
	template <typename TreeWrapper, typename RunFunc>
	static void run(commandLine &P, RunFunc test_func)
	{
		char *input_file_path = P.getOptionValue("-p");
		int K = P.getOptionIntValue("-k", 100);
		int dims = P.getOptionIntValue("-d", 3);
		size_t N = P.getOptionLongValue("-n", -1);
		int tag = P.getOptionIntValue("-t", 1);
		int rounds = P.getOptionIntValue("-r", 3);
		int query_type = P.getOptionIntValue("-q", 0);
		int read_insert_file = P.getOptionIntValue("-i", 1);
		char *insert_file_path_cml = P.getOptionValue("-I");
		int summary = P.getOptionIntValue("-s", 0);
		int tree_type = P.getOptionIntValue("-T", 0);
		int split_type = P.getOptionIntValue("-l", 0);

		using Point = typename TreeWrapper::Point;
		using points_type = parlay::sequence<Point>;
		constexpr auto num_dims = Point::get_dim();

		print_tree_param<TreeWrapper>();

		std::string name, insert_file_path = "";
		points_type wp, wi;

		if (input_file_path != NULL) { // NOTE: read main points_type
			name = std::string(input_file_path);
			name = name.substr(name.rfind('/') + 1);
			std::cout << name << " ";
			auto [n, d] =
				read_points<Point>(input_file_path, wp, 0);
			N = n;
			assert(d == num_dims);
		}

		if (read_insert_file == 1) { // NOTE: read points to be inserted
			if (insert_file_path_cml !=
			    NULL) { // given in commadnline
				insert_file_path =
					std::string(insert_file_path_cml);
				// std::cout << insert_file_path << std::endl;
			} else { // determine the name otherwise
				/*
				 * The sibling is derived from the stem, so a
				 * file not named <number>.in used to throw out
				 * of stoi, and a missing sibling used to crash
				 * after an unrelated "unable to open" message.
				 */
				std::string const stem =
					name.substr(0, name.find_first_of('.'));
				if (stem.empty() ||
				    stem.find_first_not_of("0123456789") !=
					    std::string::npos) {
					std::cerr
						<< "\n-i 1 derives the insert "
						   "file from the input's "
						   "name, "
						   "which must be <number>.in. "
						   "Pass -I <file> to name it, "
						   "or -i 0 for no insert "
						   "file.\n";
					return;
				}
				int id = std::stoi(stem);
#ifdef CCP
				id = (id + 1) %
				     10; // WARN: MOD graph number used to test
#else
				id = (id + 1) % 3;
#endif // CCP
				if (!id)
					id++;
				auto pos = std::string(input_file_path)
						   .rfind('/') +
					   1;
				insert_file_path = std::string(input_file_path)
							   .substr(0, pos) +
						   std::to_string(id) + ".in";
			}
			if (!std::filesystem::exists(insert_file_path)) {
				std::cerr << "\ninsert file "
					  << insert_file_path
					  << " does not exist. Pass -I <file> "
					     "or -i 0.\n";
				return;
			}
			auto [n, d] = read_points<Point>(
				insert_file_path.c_str(), wi, N);
			assert(d == num_dims);
		}

		// apply the test function
		test_func.template operator()<TreeWrapper, Point>(
			num_dims, wp, wi, N, K, rounds, insert_file_path, tag,
			query_type, summary);
	};

	// NOTE: For kd tree and orth tree
	template <typename RunFunc>
	static void apply_orthogonal(int const tree_type, int const dim,
				     int const split_type, commandLine &params,
				     RunFunc test_func)
	{
		auto build_tree_type = [&]<typename Point,
					   typename SplitRule>() {
			using base_type = psi::point_traits<Point>;
			if (tree_type == 0) {
				// run<kd_tree_wrapper<Point, SplitRule,
				// psi::box_leaf_aug<base_type>,
				// psi::box_interior_aug<base_type>>>(params,
				// test_func);
				run<kd_tree_wrapper<
					Point, SplitRule,
					psi::box_leaf_aug<base_type>,
					psi::box_interior_aug<base_type>>>(
					params, test_func);
			} else if (tree_type == 1) {
				/* orth_tree needs spatial_median; see the
				 * static_assert in orth_tree.h. Asking for
				 * object_median used to segfault during build.
				 */
				if constexpr (
					psi::is_spatial_median_split<
						typename SplitRule::
							partition_rule_type>) {
					run<orth_tree_wrapper<
						Point, SplitRule,
						psi::box_leaf_aug<base_type>,
						psi::box_interior_aug<
							base_type>>>(params,
								     test_func);
				} else {
					std::cout << "orth_tree needs a "
						     "spatial_median split "
						     "rule\n";
				}
			} else if (tree_type == 4) { // NOTE: for boost
				run<kd_tree_wrapper<
					Point, SplitRule,
					psi::box_leaf_aug<base_type>,
					psi::box_interior_aug<base_type>>>(
					params, test_func);
			}
		};

		// NOTE: pick the split rule
		// The lsb is the dim rule and the msb is the divide rule
		auto run_with_split_type = [&]<typename Point>() {
			if (!(split_type & (1 << 0)) &&
			    !(split_type & (1 << 1))) {
				// NOTE: 0 -> max_stretch + object_mid
				build_tree_type.template operator()<
					Point,
					psi::orthogonal_split_rule<
						psi::max_stretch_dim<Point>,
						psi::object_median<Point>>>();
			} else if ((split_type & (1 << 0)) &&
				   !(split_type & (1 << 1))) {
				// NOTE: 1 -> rotate_dim + object_mid
				build_tree_type.template operator()<
					Point,
					psi::orthogonal_split_rule<
						psi::rotate_dim<Point>,
						psi::object_median<Point>>>();
			} else if (!(split_type & (1 << 0)) &&
				   (split_type & (1 << 1))) {
				// NOTE: 2 -> max_stretch + spatial_median
				build_tree_type.template operator()<
					Point,
					psi::orthogonal_split_rule<
						psi::max_stretch_dim<Point>,
						psi::spatial_median<Point>>>();
			} else if ((split_type & (1 << 0)) &&
				   (split_type & (1 << 1))) {
				// NOTE: 3 -> rotate + spatial_median
				build_tree_type.template operator()<
					Point,
					psi::orthogonal_split_rule<
						psi::rotate_dim<Point>,
						psi::spatial_median<Point>>>();
			} else {
				std::cout << "Unsupported split type: "
					  << split_type << std::endl;
			}
		};

		if (dim == 2) {
			// run_with_split_type.template
			// operator()<basic_point<coord_type, 2>>();
			run_with_split_type.template
			operator()<aug_point<coord_type, 2, aug_id>>();
		} else if (dim == 3) {
			run_with_split_type.template
			operator()<aug_point<coord_type, 3, aug_id>>();
		}
	}

	template <typename RunFunc>
	static void
	apply_spatial_filling_curve(int const tree_type, int const dim,
				    int const split_type, commandLine &params,
				    RunFunc test_func)
	{
		auto build_tree_type =
			[&]<typename Point, typename SplitRule>() {
				if (tree_type == 0) {
					// run.template
					// operator()<kd_tree_wrapper<Point,
					// SplitRule>>();
				} else if (tree_type == 1) {
					// run.template
					// operator()<orth_tree_wrapper<Point,
					// SplitRule>>();
				} else if (tree_type == 2) {
					run<p_tree_wrapper<Point, SplitRule>>(
						params, test_func);
				} else {
					std::cout << "Unsupported tree type: "
						  << tree_type << std::endl;
				}
			};

		// NOTE: pick the split rule
		auto run_with_split_type = [&]<typename Point>() {
			if (split_type & (1 << 0)) {
				build_tree_type.template operator()<
					Point, psi::spatial_filling_curve<
						       hilbert_curve<Point>>>();
			} else if (split_type & (1 << 1)) {
				build_tree_type.template operator()<
					Point, psi::spatial_filling_curve<
						       morton_curve<Point>>>();
			}
		};

		if (dim == 2) {
			run_with_split_type.template
			operator()<aug_point<coord_type, 2, aug_id_code>>();
		} else if (dim == 3) {
			run_with_split_type.template
			operator()<aug_point<coord_type, 3, aug_id_code>>();
		}
	}

	template <typename RunFunc>
	static void apply_baselines(int const tree_type, int const dim,
				    int const split_type, commandLine &params,
				    RunFunc test_func)
	{
		auto build_tree_type = [&]<typename Point,
					   typename SplitRule>() {
			if (tree_type == 0) {
				// run.template
				// operator()<kd_tree_wrapper<Point,
				// SplitRule>>();
			} else if (tree_type == 1) {
				// run.template
				// operator()<orth_tree_wrapper<Point,
				// SplitRule>>();
			} else if (tree_type == 2) {
				// run<p_tree_wrapper<Point, SplitRule>>(params,
				// test_func);
			} else if (tree_type == 3) {
				run<cpam_raw_wrapper<Point, SplitRule>>(
					params, test_func);
			} else if (tree_type == 4) {
				; // for boost
			} else if (tree_type == 5) {
				if (dim == 2) {
					run<zd_tree_wrapper<
						typename ZD::geobase::Point,
						SplitRule>>(params, test_func);
				} else if (dim == 3) {
					run<zd_tree_3d_wrapper<
						typename ZD3D::geobase::Point,
						SplitRule>>(params, test_func);
				} else {
					std::cout << "Zdtree only supports 2D "
						     "and 3D data."
						  << std::endl;
				}
			} else {
				std::cout << "Unsupported tree type: "
					  << tree_type << std::endl;
			}
		};

		// NOTE: pick the split rule
		auto run_with_split_type = [&]<typename Point>() {
			if (split_type & (1 << 0)) {
				build_tree_type.template operator()<
					Point, psi::spatial_filling_curve<
						       hilbert_curve<Point>>>();
			} else if (split_type & (1 << 1)) {
				build_tree_type.template operator()<
					Point, psi::spatial_filling_curve<
						       morton_curve<Point>>>();
			}
		};

		if (dim == 2) {
			run_with_split_type.template
			operator()<aug_point<coord_type, 2, aug_id_code>>();
		} else if (dim == 3) {
			if (tree_type == 4) {
				;
			} else {
				run_with_split_type.template
				operator()<aug_point<coord_type, 3,
						     aug_id_code>>();
			}
		}
	}
};
