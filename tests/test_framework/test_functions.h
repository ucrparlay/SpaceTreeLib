#pragma once

#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>

#include "../common/time_loop.h"
#include "parlay/monoid.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/slice.h"
#include "psi/dependence/type_trait.h"
#include "test_config.h"
#include "test_helpers.h"

using namespace psi;

template <typename Point, typename Tree, bool kTestTime = true, int kPrint = 1>
void BuildTree(parlay::sequence<Point> const& WP, int const& rounds, Tree& pkd,
               int remaining_frac = 1) {
  using Points = typename Tree::Points;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;

  double loop_late = rounds > 1 ? 0.01 : -100;
  size_t n = WP.size();
  Points wp = Points::uninitialized(n);

  if constexpr (kTestTime) {
    pkd.DeleteTree();

    double aveBuild = time_loop(
        rounds, loop_late, [&]() { parlay::copy(WP.cut(0, n), wp.cut(0, n)); },
        [&]() { pkd.Build(wp.cut(0, n)); }, [&]() { pkd.DeleteTree(); });

    parlay::copy(WP.cut(0, n / remaining_frac), wp.cut(0, n / remaining_frac));
    pkd.Build(wp.cut(0, n / remaining_frac));

    if constexpr (kPrint == 0) {
      std::cout << aveBuild << " " << std::flush;
    } else if constexpr (kPrint == 1) {
      std::cout << aveBuild << " " << std::flush;
      if constexpr (IsKdTree<Tree> || IsOrthTree<Tree>) {
        auto deep = pkd.template GetAveTreeHeight<Leaf, Interior>();
        std::cout << deep << " " << std::flush;
      } else {
        std::cout << "-1"
                  << " " << std::flush;
      }
    } else if constexpr (kPrint == 2) {
      size_t max_deep = 0;
      std::cout << aveBuild << " ";
      if constexpr (IsKdTree<Tree> || IsOrthTree<Tree>) {
        std::cout << pkd.template GetMaxTreeDepth<Leaf, Interior>(pkd.GetRoot(),
                                                                  max_deep)
                  << " " << pkd.template GetAveTreeHeight<Leaf, Interior>()
                  << " " << std::flush;
      } else {
        std::cout << "-1 -1"
                  << " " << std::flush;
      }
    } else if constexpr (kPrint == 3) {  // for incre insert directly
      puts("# Insert");
      std::cout << "## " << 1 << std::endl;
      std::cout << "median: (1, " << aveBuild << ")-> min: (1, " << aveBuild
                << ")-> max: (1, " << aveBuild << ")-> tot: " << aveBuild
                << "-> avg: " << aveBuild << std::endl;
    }

  } else {
    // NOTE: always return a built tree
    pkd.DeleteTree();
    parlay::copy(WP.cut(0, n / remaining_frac), wp.cut(0, n / remaining_frac));
    pkd.Build(wp.cut(0, n / remaining_frac));
  }

  return;
}

template <typename Point, typename Tree, bool kTestTime = true,
          bool kSerial = false>
void BatchInsert(Tree& pkd, parlay::sequence<Point> const& WP,
                 parlay::sequence<Point> const& WI, int const& rounds,
                 double ratio = 1.0) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  using Geo = Tree::Geo;
  Points wp = Points::uninitialized(WP.size());
  Points wi = Points::uninitialized(WI.size());
  double loop_late = rounds > 1 ? 0.01 : -100;

  // NOTE: build the tree by type
  auto build_tree_by_type = [&]() {
    if constexpr (psi::IsKdTree<Tree> || psi::IsPTree<Tree>) {
      parlay::copy(WP, wp), parlay::copy(WI, wi);
      pkd.Build(parlay::make_slice(wp));
    } else if constexpr (psi::IsOrthTree<Tree>) {
      parlay::copy(WP, wp), parlay::copy(WI, wi);
      auto box1 = Geo::GetBox(parlay::make_slice(wp));
      auto box2 =
          Geo::GetBox(wi.cut(0, static_cast<size_t>(wi.size() * ratio)));
      Box box = Geo::GetBox(box1, box2);
      pkd.Build(parlay::make_slice(wp), box);
    } else {
      parlay::copy(WP, wp), parlay::copy(WI, wi);
      pkd.Build(parlay::make_slice(wp));
    }
  };

  if constexpr (kTestTime) {  // NOTE: clean and measure time
    pkd.DeleteTree();
    double aveInsert = time_loop(
        rounds, loop_late, [&]() { build_tree_by_type(); },
        [&]() {
          pkd.BatchInsert(wi.cut(0, static_cast<size_t>(wi.size() * ratio)));
        },
        [&]() { pkd.DeleteTree(); });
    std::cout << aveInsert << " " << std::flush;
    BuildTree<Point, Tree, false>(WP, rounds, pkd);
  } else {  // NOTE: insert the points from previous tree
    pkd.DeleteTree();
    build_tree_by_type();
    parlay::copy(WI, wi);
    pkd.BatchInsert(wi.cut(0, static_cast<size_t>(wi.size() * ratio)));
    std::cout << "finish" << std::endl;
  }

  return;
}

template <typename Point, typename Tree, bool kTestTime = true,
          bool kSerial = false>
void BatchDelete(Tree& pkd, parlay::sequence<Point> const& WP,
                 parlay::sequence<Point> const& WI, int const& rounds,
                 double ratio = 1.0) {
  using Points = typename Tree::Points;
  Points wp = Points::uninitialized(WP.size());
  Points wi = Points::uninitialized(WP.size());
  size_t batchSize = static_cast<size_t>(WP.size() * ratio);
  double loop_late = rounds > 1 ? 0.01 : -100;

  if constexpr (kTestTime) {
    pkd.DeleteTree();
    double aveDelete = time_loop(
        rounds, loop_late,
        [&]() {
          BuildTree<Point, Tree, false>(WP, rounds, pkd);
          parlay::copy(WI, wi);
        },
        [&]() { pkd.BatchDelete(wi.cut(0, batchSize)); },
        [&]() { pkd.DeleteTree(); });
    std::cout << aveDelete << " " << std::flush;
    BuildTree<Point, Tree, false>(WP, rounds, pkd);
  } else {
    parlay::copy(WI, wi);
    pkd.BatchDelete(wi.cut(0, batchSize));
  }

  return;
}

template <typename Point, typename Tree, bool kTestTime = true>
void BatchDiff(Tree& pkd, parlay::sequence<Point> const& WP, int const& rounds,
               double const total_ratio = 1.0,
               double const overlap_ratio = 0.5) {
  using Points = typename Tree::Points;
  size_t total_batch_size = static_cast<size_t>(WP.size() * total_ratio);
  size_t overlap_size = static_cast<size_t>(total_batch_size * overlap_ratio);

  Points WI = parlay::tabulate(total_batch_size, [&](size_t i) -> Point {
    if (i < overlap_size) {
      return WP[i];
    } else {
      Point p = WP[i];
      if constexpr (IsAugPoint<Point>) {
        p.aug.id = -1 * p.aug.id - 1;
      } else {
        std::transform(p.pnt.begin(), p.pnt.end(), p.pnt.begin(),
                       std::negate<Coord>());
      }
      return p;
    }
  });

  Points wp = Points::uninitialized(WP.size());
  Points wi = Points::uninitialized(WI.size());

  if constexpr (kTestTime) {
    pkd.DeleteTree();
    double aveDelete = time_loop(
        rounds, 0.01,
        [&]() {
          BuildTree<Point, Tree, false>(WP, rounds, pkd);
          parlay::copy(WI, wi);
        },
        [&]() { pkd.BatchDiff(wi.cut(0, total_batch_size)); },
        [&]() { pkd.DeleteTree(); });
    std::cout << aveDelete << " " << std::flush;
    BuildTree<Point, Tree, false>(WP, rounds, pkd);
  } else {
    parlay::copy(WI, wi);
    pkd.BatchDiff(wi.cut(0, total_batch_size));
  }

  return;
}

template <typename Point, typename Tree, bool kInsert>
void BatchInsertByStep(Tree& pkd, parlay::sequence<Point> const& WP,
                       int const rounds, double const insert_ratio,
                       int const remain_divide_ratio = 2) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  using Geo = Tree::Geo;
  Points wp = Points::uninitialized(WP.size());
  size_t n = wp.size();
  size_t step = static_cast<size_t>(insert_ratio * n);
  size_t slice_num = n / step;
  parlay::sequence<parlay::sequence<double>> time_table(
      rounds + 2, parlay::sequence<double>(slice_num, 0.0));
  parlay::sequence<StepUpdateLogger> log_time(slice_num);
  int id_cnt = 0;
  for (auto& i : log_time) {
    i = {id_cnt++, 0.0};
  }
  size_t round_cnt = 0;

  pkd.DeleteTree();

  // NOTE: build the tree by type
  auto prepare_build = [&]() {
    pkd.DeleteTree();
    if constexpr (psi::IsKdTree<Tree> || psi::IsPTree<Tree>) {
      parlay::copy(WP, wp);
    } else if constexpr (psi::IsOrthTree<Tree>) {
      parlay::copy(WP, wp);
      auto box = Geo::GetBox(wp.cut(0, n));
      pkd.SetBoundingBox(box);
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
      pkd.BatchInsert(parlay::make_slice(wp.begin() + l, wp.begin() + r));
      l = r;
      time_table[round_cnt][cnt++] += t.next_time();
    }
    round_cnt++;
  };

  double loop_late = rounds > 1 ? 0.01 : -100;
  double ave_time = time_loop(
      rounds, loop_late, [&]() { prepare_build(); },
      [&]() { incre_build(slice_num + 1); }, [&]() { pkd.DeleteTree(); });

  // begin count the time
  if (rounds != 1 && round_cnt - 1 != rounds) {
    throw std::runtime_error("rounds not match!");
  }
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
            [](auto const& a, auto const& b) { return a.t < b.t; });

  // Calculate average time from log_time
  double total_time = 0.0;
  for (auto const& log : log_time) {
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

  return;
}

template <typename Point, typename Tree, bool kInsert>
void BatchDeleteByStep(Tree& pkd, parlay::sequence<Point> const& WP,
                       int const rounds, double const insert_ratio,
                       size_t const remain_divide_ratio = 2) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  using Geo = Tree::Geo;
  Points wp = Points::uninitialized(WP.size());
  size_t n = wp.size();
  size_t step = static_cast<size_t>(insert_ratio * n);
  size_t slice_num = n / step;
  parlay::sequence<parlay::sequence<double>> time_table(
      rounds + 2, parlay::sequence<double>(slice_num, 0.0));
  parlay::sequence<StepUpdateLogger> log_time(slice_num);
  int id_cnt = 0;
  for (auto& i : log_time) {
    i = {id_cnt++, 0.0};
  }
  size_t round_cnt = 0;

  pkd.DeleteTree();

  // NOTE: build the tree by type
  auto build_tree_by_type = [&]() {
    parlay::copy(WP, wp);

    if constexpr (psi::IsKdTree<Tree> || psi::IsPTree<Tree>) {
      pkd.Build(parlay::make_slice(wp));
    } else if constexpr (psi::IsOrthTree<Tree>) {
      auto box = Geo::GetBox(wp.cut(0, n));
      pkd.Build(parlay::make_slice(wp), box);
    } else {
      pkd.Build(parlay::make_slice(wp));
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
      pkd.BatchDelete(
          parlay::make_slice(wp_slice.begin() + l, wp_slice.begin() + r));
      l = r;
      time_table[round_cnt][cnt++] += t.next_time();
    }
    round_cnt++;
  };

  double loop_late = rounds > 1 ? 0.01 : -100;
  double ave_time = time_loop(
      rounds, loop_late, [&]() { build_tree_by_type(); },
      [&]() { incre_delete(slice_num, wp.cut(0, slice_num * step)); },
      [&]() { pkd.DeleteTree(); });

  if (rounds != 1 && round_cnt - 1 != rounds) {
    throw std::runtime_error("rounds not match!");
  }
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
  std::sort(log_time.begin(), log_time.end(),
            [](auto const& a, auto const& b) { return a.t < b.t; });
  // Calculate average time from log_time
  double total_time = 0.0;
  for (auto const& log : log_time) {
    total_time += log.t;
  }
  double average_time = total_time / log_time.size();

  std::cout << "median: " << log_time[slice_num / 2]
            << "-> min: " << *log_time.begin()
            << "-> max: " << *log_time.rbegin() << "-> tot: " << ave_time
            << "-> avg: " << average_time << std::endl;

  // WARN: restore status
  pkd.DeleteTree();
  build_tree_by_type();
  incre_delete(slice_num / remain_divide_ratio,
               wp.cut(slice_num / remain_divide_ratio * step, n));

  return;
}

template <typename Point, typename Tree, bool printHeight = 0,
          bool printVisNode = 1>
void QueryKNN([[maybe_unused]] uint_fast8_t const& Dim,
              parlay::sequence<Point> const& WP, int const& rounds, Tree& pkd,
              Typename* kdknn, int const K, bool const flattenTreeTag) {
  using Points = typename Tree::Points;
  using Coord = typename Point::Coord;
  using DisType = typename Point::DisType;
  using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  size_t n = WP.size();
  double loopLate = rounds > 1 ? 0.01 : -0.1;
  auto* KDParallelRoot = pkd.GetRoot();

  Points wp = WP;

  parlay::sequence<nn_pair> Out(
      K * n, nn_pair(std::ref(wp[0]), static_cast<DisType>(0)));
  parlay::sequence<kBoundedQueue<Point, nn_pair>> bq =
      parlay::sequence<kBoundedQueue<Point, nn_pair>>::uninitialized(n);
  parlay::parallel_for(
      0, n, [&](size_t i) { bq[i].resize(Out.cut(i * K, i * K + K)); });
  parlay::sequence<size_t> vis_leaf(n), vis_inter(n), gen_box(n), check_box(n),
      skip_box(n);

  double aveQuery = time_loop(
      rounds, loopLate,
      [&]() { parlay::parallel_for(0, n, [&](size_t i) { bq[i].reset(); }); },
      [&]() {
        parlay::parallel_for(0, n, [&](size_t i) {
          auto [vis_leaf_num, vis_inter_num, gen_box_num, check_box_num,
                skip_box_num] = pkd.KNN(KDParallelRoot, wp[i], bq[i]);
          kdknn[i] = bq[i].top().second;
          vis_leaf[i] = vis_leaf_num;
          vis_inter[i] = vis_inter_num;
          gen_box[i] = gen_box_num;
          check_box[i] = check_box_num;
          skip_box[i] = skip_box_num;
        });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;
  if (printHeight) {
    size_t max_deep = 0;
    std::cout << pkd.template GetMaxTreeDepth<Leaf, Interior>(pkd.GetRoot(),
                                                              max_deep)
              << " " << pkd.template GetAveTreeHeight<Leaf, Interior>() << " "
              << std::flush;
  }
  if (printVisNode) {
    std::cout << static_cast<double>(parlay::reduce(vis_leaf.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(vis_inter.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(gen_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(check_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(skip_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
  }

  return;
}

// NOTE: run range count and check the correct
template <typename Point, typename Tree>
void RangeCount(parlay::sequence<Point> const& wp, Tree& pkd, int const& rounds,
                int rec_num, int rec_type, int const DIM) {
  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, false>(rec_num, rec_type, wp, DIM);
  parlay::sequence<size_t> kdknn(rec_num, 0);

  double aveCount = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeCount(query_box_seq[i].first);
          kdknn[i] = size;
        });
      },
      [&]() {});

  // NOTE: verify the solutions
  std::cout << "check range count: " << rec_num << " " << rec_type << std::endl;
  for (int i = 0; i < rec_num; i++) {
    if (kdknn[i] != query_box_seq[i].second) {
      std::cout << kdknn[i] << " " << query_box_seq[i].second << " "
                << query_box_seq[i].first.first << " "
                << query_box_seq[i].first.second << std::endl;
    }
    assert(std::cmp_equal(kdknn[i], query_box_seq[i].second));
  }

  std::cout << aveCount << "\n" << std::flush;

  return;
}

template <typename Point, typename Tree>
void rangeCountRadius(parlay::sequence<Point> const& wp, Tree& pkd,
                      Typename* kdknn, int const& rounds, int const& queryNum) {
  using Box = typename Tree::Box;
  using circle = typename Tree::circle;

  int n = wp.size();

  double aveCount = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, queryNum, [&](size_t i) {
          Box query_box_seq = pkd.get_box(
              Box(wp[i], wp[i]), Box(wp[(i + n / 2) % n], wp[(i + n / 2) % n]));
          auto d =
              Tree::p2p_distance(wp[i], wp[(i + n / 2) % n], wp[i].GetDim());
          d = static_cast<Coord>(std::sqrt(d));
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
void RangeQuery(parlay::sequence<Point> const& wp, Tree& pkd, int const& rounds,
                int const rec_num, int const rec_type, int const DIM) {
  using Points = typename Tree::Points;

  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, true>(rec_num, rec_type, wp, DIM);
  auto [offset, tot_size] = parlay::scan(
      parlay::delayed_tabulate(
          rec_num,
          [&](size_t i) -> size_t { return query_box_seq[i].second.size(); }),
      parlay::addm<size_t>());
  offset.push_back(tot_size);
  Points Out(tot_size);
  parlay::sequence<size_t> kdknn(rec_num, 0);
  parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
      gen_box(rec_num), full_box(rec_num), skip_box(rec_num);

  double aveQuery = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeQuery(
              query_box_seq[i].first, Out.cut(offset[i], offset[i + 1]));
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
      std::cout << kdknn[i] << " " << query_box_seq[i].second.size() << " "
                << query_box_seq[i].first.first << query_box_seq[i].first.second
                << std::endl;
      std::cout << vis_leaf[i] << " " << gen_box[i] << " " << full_box[i] << " "
                << skip_box[i] << std::endl;
    }
    assert(std::cmp_equal(kdknn[i], query_box_seq[i].second.size()));
    parlay::sort_inplace(
        Out.cut(offset[i], offset[i + 1]),
        [&](auto const& a, auto const& b) { return a.aug.id < b.aug.id; });
    parlay::sort_inplace(
        query_box_seq[i].second,
        [&](auto const& a, auto const& b) { return a.aug.id < b.aug.id; });
    for (size_t j = 0; j < query_box_seq[i].second.size(); j++) {
      if (Out[offset[i] + j] != query_box_seq[i].second.at(j)) {
        std::cout << "wrong" << query_box_seq[i].first.first
                  << query_box_seq[i].first.second << std::endl;
        std::cout << Out[offset[i] + j] << " " << query_box_seq[i].second.at(j)
                  << std::endl;
      }

      if constexpr (IsKdTree<Tree> ||
                    IsOrthTree<Tree>) {  // TODO: fix this by enable kdtree
                                         // handling duplicates by id
        assert(Out[offset[i] + j] == query_box_seq[i].second.at(j));
      } else if constexpr (IsPTree<Tree>) {
        assert(Out[offset[i] + j] == query_box_seq[i].second.at(j));
      }
    }
  }

  std::cout << aveQuery << "\n" << std::flush;
  return;
}

//* test range count for fix rectangle
template <typename Point, typename Tree>
void RangeCountFix(Tree& pkd, Typename* kdknn, int const& rounds, int rec_type,
                   int rec_num, int DIM, auto const& query_box_seq,
                   auto max_size) {
  parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
      gen_box(rec_num), full_box(rec_num), skip_box(rec_num);

  double aveCount = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeCount(query_box_seq[i].first);

          kdknn[i] = size;
          vis_leaf[i] = logger.vis_leaf_num;
          vis_inter[i] = logger.vis_interior_num;
          gen_box[i] = logger.generate_box_num;
          full_box[i] = logger.full_box_num;
          skip_box[i] = logger.skip_box_num;
        });
      },
      [&]() {});

  std::cout << aveCount << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_leaf.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_inter.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(gen_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(full_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(skip_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;

  return;
}

//* test range query for fix rectangle
template <typename Point, typename Tree>
void RangeQueryFix(Tree& pkd, Typename* kdknn, int const& rounds,
                   parlay::sequence<Point>& Out, int rec_type, int rec_num,
                   int DIM, auto const& query_box_seq, auto max_size) {
  parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
      gen_box(rec_num), full_box(rec_num), skip_box(rec_num);
  auto [offset, tot_size] = parlay::scan(
      parlay::delayed_tabulate(
          rec_num, [&](size_t i) -> size_t { return query_box_seq[i].second; }),
      parlay::addm<size_t>());
  offset.push_back(tot_size);
  Out.resize(tot_size);

  double aveQuery = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeQuery(
              query_box_seq[i].first, Out.cut(offset[i], offset[i + 1]));

          kdknn[i] = size;
          vis_leaf[i] = logger.vis_leaf_num;
          vis_inter[i] = logger.vis_interior_num;
          gen_box[i] = logger.generate_box_num;
          full_box[i] = logger.full_box_num;
          skip_box[i] = logger.skip_box_num;

          assert(std::cmp_equal(kdknn[i], query_box_seq[i].second));
        });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_leaf.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_inter.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(gen_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(full_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(skip_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  return;
}

template <typename Point, typename Tree>
void RangeQuerySerialWithLog(Tree& pkd, Typename* kdknn, int const& rounds,
                             parlay::sequence<Point>& Out, int rec_type,
                             int rec_num, int DIM, auto const& query_box_seq,
                             auto max_size) {
  auto [offset, tot_size] = parlay::scan(
      parlay::delayed_tabulate(
          rec_num, [&](size_t i) -> size_t { return query_box_seq[i].second; }),
      parlay::addm<size_t>());
  offset.push_back(tot_size);
  Out.resize(tot_size);

  for (int i = 0; i < rec_num; i++) {
    parlay::internal::timer t;
    t.reset(), t.start();
    auto [size, logger] = pkd.RangeQuery(query_box_seq[i].first,
                                         Out.cut(offset[i], offset[i + 1]));
    t.stop();
    std::cout << rec_type << " " << query_box_seq[i].second << " "
              << std::scientific << t.total_time() << std::endl;
  }

  return;
}
