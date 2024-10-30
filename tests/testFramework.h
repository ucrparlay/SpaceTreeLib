#pragma once

#include <algorithm>
#include <boost/random/linear_feedback_shift.hpp>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ios>

#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "common/time_loop.h"
#include "cpdd/base_tree.h"
#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "cpdd/orth_tree.h"
#include "parlay/internal/group_by.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/slice.h"

using Coord = long;
// using Coord = double;
using Typename = Coord;
using namespace cpdd;

// NOTE: KNN size
static constexpr double batchQueryRatio = 0.01;
static constexpr size_t batchQueryOsmSize = 10000000;
// NOTE: rectangle numbers
static constexpr int rangeQueryNum = 50000;
static constexpr int singleQueryLogRepeatNum = 100;

// NOTE: rectangle numbers for inba ratio
static constexpr int rangeQueryNumInbaRatio = 50000;
// NOTE: insert batch ratio for inba ratio
static constexpr double insertBatchInbaRatio = 0.001;
// NOTE: knn batch ratio for inba ratio
static constexpr double knnBatchInbaRatio = 0.1;

// NOTE: Insert Ratio when summary
static constexpr double batchInsertRatio = 0.01;
// NOTE: rectange type used in summary
static constexpr int summaryRangeQueryType = 2;
// NOTE: range query num in summary
static constexpr int summaryRangeQueryNum = 10000;

// NOTE: helper for delete type
enum DeleteType { kBatchDelete, kBatchDiff };

// * [a,b)
inline size_t get_random_index(size_t a, size_t b, [[maybe_unused]] int seed) {
  return size_t((rand() % (b - a)) + a);
  // return size_t( ( parlay::hash64( static_cast<uint64_t>( seed ) ) % ( b
  // - a
  // ) ) + a );
}

template <typename Point, typename Tree, bool SavePoint>
size_t recurse_box(parlay::slice<Point*, Point*> In, auto& box_seq, int DIM,
                   std::pair<size_t, size_t> range, int& idx, int rec_num,
                   int type) {
  // using Tree = Tree;
  // using Box = typename Tree::Box;

  size_t n = In.size();
  if (idx >= rec_num || n < range.first || n == 0) return 0;

  size_t mx = 0;
  bool goon = false;
  if (n <= range.second) {
    if constexpr (SavePoint) {
      box_seq[idx++] =
          std::make_pair(Tree::GetBox(In), parlay::to_sequence(In));
    } else {
      box_seq[idx++] = std::make_pair(Tree::GetBox(In), In.size());
    }

    // NOTE: handle the cose that all Points are the same then become un -
    // divideable
    // Modify the coefficient to make the rectangle size
    // distribute as uniform as possible within the range
    if ((type == 0 && n < 40 * range.first) ||
        (type == 1 && n < 10 * range.first) ||
        (type == 2 && n < 2 * range.first) ||
        parlay::all_of(In, [&](Point const& p) { return p == In[0]; })) {
      return In.size();
    } else {
      goon = true;
      mx = n;
    }
  }

  int dim = get_random_index(0, DIM, rand());
  size_t pos = get_random_index(0, n, rand());
  parlay::sequence<bool> flag(n, 0);
  parlay::parallel_for(0, n, [&](size_t i) {
    if (cpdd::Num_Comparator<Coord>::Gt(In[i].pnt[dim], In[pos].pnt[dim]))
      flag[i] = 1;
    else
      flag[i] = 0;
  });
  auto [Out, m] = parlay::internal::split_two(In, flag);

  assert(Out.size() == n);
  // LOG << dim << " " << Out[0] << Out[m] << ENDL;
  size_t l, r;
  l = recurse_box<Point, Tree, SavePoint>(Out.cut(0, m), box_seq, DIM, range,
                                          idx, rec_num, type);
  r = recurse_box<Point, Tree, SavePoint>(Out.cut(m, n), box_seq, DIM, range,
                                          idx, rec_num, type);

  if (goon) {
    return mx;
  } else {
    return std::max(l, r);
  }
}

template <typename Point, typename Tree, bool SavePoint>
auto gen_rectangles(int rec_num, int const type,
                    parlay::sequence<Point> const& WP, int DIM) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  using BoxSeq = std::conditional_t<SavePoint == false,
                                    parlay::sequence<std::pair<Box, size_t>>,
                                    parlay::sequence<std::pair<Box, Points>>>;

  size_t n = WP.size();
  std::pair<size_t, size_t> range;
  if (type == 0) {  //* small bracket
    range.first = 1;
    range.second = size_t(std::sqrt(std::sqrt(n)));
  } else if (type == 1) {  //* medium bracket
    range.first = size_t(std::sqrt(std::sqrt(n)));
    range.second = size_t(std::sqrt(n));
  } else if (type == 2) {  //* large bracket
    range.first = size_t(std::sqrt(n));

    // NOTE: special handle for large dimension datasets
    if (n == 100000000)
      range.second = n / 100 - 1;
    else if (n == 1000000000)
      range.second = n / 1000 - 1;
    else
      range.second = n - 1;
  }
  BoxSeq box_seq(rec_num);
  int cnt = 0;
  Points wp(n);

  srand(10);

  // LOG << " " << range.first << " " << range.second << ENDL;

  size_t max_size = 0;
  while (cnt < rec_num) {
    parlay::copy(WP, wp);
    auto r = recurse_box<Point, Tree, SavePoint>(
        parlay::make_slice(wp), box_seq, DIM, range, cnt, rec_num, type);
    max_size = std::max(max_size, r);
    // LOG << cnt << " " << max_size << ENDL;
  }
  // LOG << "finish generate " << ENDL;
  return std::make_pair(box_seq, max_size);
}

template <typename Point, typename Tree, int kPrint = 1>
void buildTree([[maybe_unused]] int const& Dim,
               parlay::sequence<Point> const& WP, int const& rounds,
               Tree& pkd) {
  using Points = typename Tree::Points;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;

  double loopLate = rounds > 1 ? 1.0 : -0.1;
  size_t n = WP.size();
  Points wp = Points::uninitialized(n);
  pkd.DeleteTree();
  double aveBuild = time_loop(
      rounds, loopLate, [&]() { parlay::copy(WP.cut(0, n), wp.cut(0, n)); },
      [&]() { pkd.Build(wp.cut(0, n)); }, [&]() { pkd.DeleteTree(); });

  //* return a built Tree
  parlay::copy(WP.cut(0, n), wp.cut(0, n));
  pkd.Build(wp.cut(0, n));

  if (kPrint == 1) {
    LOG << aveBuild << " " << std::flush;
    auto deep = pkd.template GetAveTreeHeight<Leaf, Interior>();
    LOG << deep << " " << std::flush;
  } else if (kPrint == 2) {
    size_t max_deep = 0;
    LOG << aveBuild << " ";
    LOG << pkd.template GetMaxTreeDepth<Leaf, Interior>(pkd.GetRoot(), max_deep)
        << " " << pkd.template GetAveTreeHeight<Leaf, Interior>() << " "
        << std::flush;
  }

  return;
}
//
// template<typename Point, int print = 1>
// void incrementalBuild(const int Dim, const parlay::sequence<Point>& WP, const
// int rounds, BaseTree<Point>& pkd,
//                       double stepRatio) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     size_t n = WP.size();
//     size_t step = n * stepRatio;
//     Points wp = Points::uninitialized(n);
//
//     pkd.delete_tree();
//     double delay = 1.0;
//     if (pkd.get_imbalance_ratio() == 1) {
//         delay = -0.1;
//     }
//
//     double aveIncreBuild = time_loop(
//         rounds, delay, [&]() { parlay::copy(WP, wp); },
//         [&]() {
//             size_t l = 0, r = 0;
//             while (l < n) {
//                 r = std::min(l + step, n);
//                 pkd.BatchInsert(wp.cut(l, r), Dim);
//                 l = r;
//             }
//         },
//         [&]() { pkd.delete_tree(); });
//
//     parlay::copy(WP, wp);
//     size_t l = 0, r = 0;
//     while (l < n) {
//         r = std::min(l + step, n);
//         pkd.BatchInsert(wp.cut(l, r), Dim);
//         l = r;
//     }
//
//     if (print == 1) {
//         auto deep = pkd.getAveTreeHeight();
//         LOG << aveIncreBuild << " " << deep << " " << std::flush;
//     } else if (print == 2) {  // NOTE: print the maxtree height and avetree
//     height
//         size_t max_deep = 0;
//         LOG << aveIncreBuild << " " << pkd.getMaxTreeDepth(pkd.get_root(),
//         max_deep) << " " << pkd.getAveTreeHeight()
//             << " " << std::flush;
//     }
//     return;
// }
//
// template<typename Point, bool print = 1>
// void incrementalDelete(const int Dim, const parlay::sequence<Point>& WP,
// const parlay::sequence<Point>& WI, int rounds,
//                        BaseTree<Point>& pkd, double stepRatio) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     size_t n = WP.size();
//     size_t step = n * stepRatio;
//     Points wp = Points::uninitialized(2 * n);
//     Points wi = Points::uninitialized(n);
//
//     pkd.delete_tree();
//
//     double aveIncreDelete = time_loop(
//         rounds, 1.0,
//         [&]() {
//             parlay::copy(WP.cut(0, n), wp.cut(0, n));
//             parlay::copy(WI.cut(0, n), wp.cut(n, 2 * n));
//             parlay::copy(WI.cut(0, n), wi.cut(0, n));
//
//             pkd.build(wp.cut(0, 2 * n), Dim);
//         },
//         [&]() {
//             size_t l = 0, r = 0;
//             while (l < n) {
//                 r = std::min(l + step, n);
//                 pkd.batchDelete(wi.cut(l, r), Dim);
//                 l = r;
//             }
//         },
//         [&]() { pkd.delete_tree(); });
//
//     parlay::copy(WP.cut(0, n), wp.cut(0, n));
//     parlay::copy(WI.cut(0, n), wp.cut(n, 2 * n));
//     parlay::copy(WI.cut(0, n), wi.cut(0, n));
//     pkd.build(wp.cut(0, 2 * n), Dim);
//     size_t l = 0, r = 0;
//     while (l < n) {
//         r = std::min(l + step, n);
//         pkd.batchDelete(wi.cut(l, r), Dim);
//         l = r;
//     }
//
//     if (print) {
//         auto deep = pkd.getAveTreeHeight();
//         LOG << aveIncreDelete << " " << deep << " " << std::flush;
//     }
//     return;
// }
//
template <typename Point, typename Tree, bool serial = false>
void BatchInsert(Tree& pkd, parlay::sequence<Point> const& WP,
                 parlay::sequence<Point> const& WI,
                 [[maybe_unused]] uint_fast8_t const& DIM, int const& rounds,
                 double ratio = 1.0) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  Points wp = Points::uninitialized(WP.size());
  Points wi = Points::uninitialized(WI.size());

  pkd.DeleteTree();

  double aveInsert = time_loop(
      rounds, 1.0,
      [&]() {
        if constexpr (cpdd::IsKdTree<Tree>) {
          parlay::copy(WP, wp), parlay::copy(WI, wi);
          pkd.Build(parlay::make_slice(wp));
        } else if constexpr (cpdd::IsOrthTree<Tree>) {
          parlay::copy(WP, wp), parlay::copy(WI, wi);
          Box Box = Tree::GetBox(Tree::GetBox(parlay::make_slice(wp)),
                                 Tree::GetBox(parlay::make_slice(wi)));
          pkd.Build(parlay::make_slice(wp), Box);
        } else {
          LOG << "Not supported Tree type" << ENDL;
        }
      },
      [&]() { pkd.BatchInsert(wi.cut(0, size_t(wi.size() * ratio))); },
      [&]() { pkd.DeleteTree(); });

  if constexpr (cpdd::IsKdTree<Tree>) {
    parlay::copy(WP, wp), parlay::copy(WI, wi);
    pkd.Build(parlay::make_slice(wp));
  } else if constexpr (cpdd::IsOrthTree<Tree>) {
    parlay::copy(WP, wp), parlay::copy(WI, wi);
    Box Box = Tree::GetBox(Tree::GetBox(parlay::make_slice(wp)),
                           Tree::GetBox(wi.cut(0, size_t(wi.size() * ratio))));
    pkd.Build(parlay::make_slice(wp), Box);
  } else {
    LOG << "Not supported Tree type" << ENDL;
  }

  pkd.BatchInsert(wi.cut(0, size_t(wi.size() * ratio)));

  LOG << aveInsert << " " << std::flush;

  return;
}

template <typename Point, typename Tree, DeleteType kDeleteType = kBatchDelete>
void batchDelete(Tree& pkd, parlay::sequence<Point> const& WP,
                 parlay::sequence<Point> const& WI,
                 [[maybe_unused]] uint_fast8_t const& DIM, int const& rounds,
                 bool afterInsert = 1, double ratio = 1.0) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  Points wp = Points::uninitialized(WP.size());
  Points wi = Points::uninitialized(
      WP.size());  //! warnning need to adjust space if necessary

  pkd.DeleteTree();
  size_t batchSize = static_cast<size_t>(WP.size() * ratio);

  double aveDelete = time_loop(
      rounds, 1.0,
      [&]() {
        if (afterInsert) {  //* first insert wi then delete wi
          parlay::copy(WP, wp), parlay::copy(WI, wi);
          if constexpr (cpdd::IsKdTree<Tree>) {
            pkd.Build(parlay::make_slice(wp));
          } else if constexpr (cpdd::IsOrthTree<Tree>) {
            Box Box = Tree::GetBox(Tree::GetBox(parlay::make_slice(wp)),
                                   Tree::GetBox(wi.cut(0, batchSize)));
            pkd.Build(parlay::make_slice(wp), Box);
          }
          pkd.BatchInsert(wi.cut(0, batchSize));
          parlay::copy(WP, wp), parlay::copy(WI, wi);
        } else {  //* only build wp and then delete from wp
          parlay::copy(WP, wp), parlay::copy(WP, wi);
          pkd.Build(parlay::make_slice(wp));
          // parlay::copy(WP, wp), parlay::copy(WP, wi);
        }
      },
      [&]() {
        if constexpr (kDeleteType == kBatchDelete) {
          pkd.BatchDelete(wi.cut(0, batchSize));
        } else if constexpr (kDeleteType == kBatchDiff) {
          pkd.BatchDiff(wi.cut(0, batchSize));
        }
      },
      [&]() { pkd.DeleteTree(); });

  if (afterInsert) {
    //* first insert wi then delete wi
    parlay::copy(WP, wp), parlay::copy(WI, wi);
    if constexpr (cpdd::IsKdTree<Tree>) {
      pkd.Build(parlay::make_slice(wp));
    } else if constexpr (cpdd::IsOrthTree<Tree>) {
      Box Box = Tree::GetBox(Tree::GetBox(parlay::make_slice(wp)),
                             Tree::GetBox(wi.cut(0, batchSize)));
      pkd.Build(parlay::make_slice(wp), Box);
    }
    pkd.BatchInsert(wi.cut(0, batchSize));
    parlay::copy(WP, wp), parlay::copy(WI, wi);
    if constexpr (kDeleteType == kBatchDelete) {
      pkd.BatchDelete(wi.cut(0, batchSize));
    } else if constexpr (kDeleteType == kBatchDiff) {
      pkd.BatchDiff(wi.cut(0, batchSize));
    }
  } else {
    parlay::copy(WP, wp);
    pkd.Build(parlay::make_slice(wp));
  }

  std::cout << aveDelete << " " << std::flush;

  return;
}

// template<typename Point, bool insert>
// void batchUpdateByStep(BaseTree<Point>& pkd, const parlay::sequence<Point>&
// WP, const parlay::sequence<Point>& WI,
//                        const uint_fast8_t& DIM, const int& rounds, double
//                        ratio = 1.0, double max_ratio = 1e-2) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     Points wp = Points::uninitialized(WP.size());
//     Points wi = Points::uninitialized(WI.size());
//     size_t n = static_cast<size_t>(max_ratio * wi.size());
//
//     pkd.delete_tree();
//
//     double aveInsert = time_loop(
//         rounds, 1.0,
//         [&]() {
//             parlay::copy(WP, wp), parlay::copy(WI, wi);
//             pkd.build(parlay::make_slice(wp), DIM);
//         },
//         [&]() {
//             size_t l = 0, r = 0;
//             size_t step = static_cast<size_t>(wi.size() * ratio);
//             while (l < n) {
//                 r = std::min(l + step, n);
//                 // LOG << l << ' ' << r << ENDL;
//                 if (insert) {
//                     pkd.BatchInsert(parlay::make_slice(wi.begin() + l,
//                     wi.begin() + r), DIM);
//                 } else {
//                     pkd.batchDelete(parlay::make_slice(wi.begin() + l,
//                     wi.begin() + r), DIM);
//                 }
//                 l = r;
//             }
//             // LOG << l << ENDL;
//         },
//         [&]() { pkd.delete_tree(); });
//
//     // WARN: not reset status
//
//     // parlay::copy(WP, wp), parlay::copy(WI, wi);
//     // pkd.build(parlay::make_slice(wp), DIM);
//     // size_t l = 0, r = 0;
//     // size_t step = wi.size() * ratio;
//     // while (l < n) {
//     //     r = std::min(l + step, n);
//     //     if (insert) {
//     //         pkd.BatchInsert(parlay::make_slice(wi.begin() + l, wi.begin()
//     + r), DIM);
//     //     } else {
//     //         pkd.batchDelete(parlay::make_slice(wi.begin() + l, wi.begin()
//     + r), DIM);
//     //     }
//     //     l = r;
//     // }
//
//     LOG << aveInsert << " " << std::flush;
//
//     return;
// }
//
template <typename Point, typename Tree, bool printHeight = 0,
          bool printVisNode = 1>
void queryKNN([[maybe_unused]] uint_fast8_t const& Dim,
              parlay::sequence<Point> const& WP, int const& rounds, Tree& pkd,
              Typename* kdknn, int const K, bool const flattenTreeTag) {
  using Points = typename Tree::Points;
  using Coord = typename Point::Coord;
  using nn_pair = std::pair<std::reference_wrapper<Point>, Coord>;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  // using nn_pair = std::pair<Point, Coord>;
  size_t n = WP.size();
  // int LEAVE_WRAP = 32;
  double loopLate = rounds > 1 ? 1.0 : -0.1;
  auto* KDParallelRoot = pkd.GetRoot();
  // auto* old = pkd.GetRoot();
  // if constexpr (cpdd::IsMultiNode<Interior>) {
  //     // LOG << "start expanding the Tree" << ENDL;
  //     KDParallelRoot =
  //         pkd.template Expand2Binary<typename Tree::KdInteriorNode,
  //         Interior>(
  //             KDParallelRoot);
  //     pkd.SetRoot(KDParallelRoot);
  //     // LOG << "finish expanding the Tree" << ENDL;
  //     // pkd.template Validate<Leaf, typename Tree::KdInteriorNode,
  //     //                       typename Tree::SplitRuleType>(Dim);
  // }

  Points wp = Points::uninitialized(n);
  parlay::copy(WP, wp);

  parlay::sequence<nn_pair> Out(K * n, nn_pair(std::ref(wp[0]), 0));
  // parlay::sequence<nn_pair> Out(K * n);
  parlay::sequence<kBoundedQueue<Point, nn_pair>> bq =
      parlay::sequence<kBoundedQueue<Point, nn_pair>>::uninitialized(n);
  parlay::parallel_for(
      0, n, [&](size_t i) { bq[i].resize(Out.cut(i * K, i * K + K)); });
  parlay::sequence<size_t> vis_nodes(n), gen_box(n), check_box(n), skip_box(n);

  double aveQuery = time_loop(
      rounds, loopLate,
      [&]() { parlay::parallel_for(0, n, [&](size_t i) { bq[i].reset(); }); },
      [&]() {
        if (!flattenTreeTag) {  // WARN: Need ensure pkd.size() == wp.size()
          pkd.Flatten(parlay::make_slice(wp));
        }
        parlay::parallel_for(0, n, [&](size_t i) {
          auto [vis_node_num, gen_box_num, check_box_num, skip_box_num] =
              pkd.KNN(KDParallelRoot, wp[i], bq[i]);
          kdknn[i] = bq[i].top().second;
          vis_nodes[i] = vis_node_num;
          gen_box[i] = gen_box_num;
          check_box[i] = check_box_num;
          skip_box[i] = skip_box_num;
        });
      },
      [&]() {});

  LOG << aveQuery << " " << std::flush;
  if (printHeight) {
    // WARN: change when using multi-node
    size_t max_deep = 0;
    LOG << pkd.template GetMaxTreeDepth<Leaf, Interior>(pkd.GetRoot(), max_deep)
        << " " << pkd.template GetAveTreeHeight<Leaf, Interior>() << " "
        << std::flush;
    // size_t max_deep = 0;
    // LOG << pkd.template GetMaxTreeDepth<Leaf,
    //                                     typename Tree::KdInteriorNode>(
    //            pkd.GetRoot(), max_deep)
    //     << " "
    //     << pkd.template GetAveTreeHeight<Leaf,
    //                                      typename Tree::KdInteriorNode>()
    //     << " " << std::flush;
  }
  if (printVisNode) {
    LOG << parlay::reduce(vis_nodes.cut(0, n)) / n << " " << std::flush;
    LOG << parlay::reduce(gen_box.cut(0, n)) / n << " " << std::flush;
    LOG << parlay::reduce(check_box.cut(0, n)) / n << " " << std::flush;
    LOG << parlay::reduce(skip_box.cut(0, n)) / n << " " << std::flush;
  }
  // pkd.SetRoot(old);

  return;
}

// NOTE: run range count and check the correct
template <typename Point, typename Tree>
void rangeCount(parlay::sequence<Point> const& wp, Tree& pkd, Typename* kdknn,
                int const& rounds, int rec_num, int rec_type, int const DIM) {
  // using Points = typename Tree::Points;
  // using Box = typename Tree::Box;

  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, false>(rec_num, rec_type, wp, DIM);

  double aveCount = time_loop(
      rounds, 1.0, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeCount(query_box_seq[i].first);
          kdknn[i] = size;
        });
      },
      [&]() {});

  // NOTE: verify the solutions
  LOG << "check range count: " << rec_num << " " << rec_type << ENDL;
  for (int i = 0; i < rec_num; i++) {
    assert(std::cmp_equal(kdknn[i], query_box_seq[i].second));
  }

  LOG << aveCount << " " << std::flush;

  return;
}

template <typename Point, typename Tree>
void rangeCountRadius(parlay::sequence<Point> const& wp, Tree& pkd,
                      Typename* kdknn, int const& rounds, int const& queryNum) {
  // using Points = typename Tree::Points;
  // using node = typename Tree::node;
  using Box = typename Tree::Box;
  using circle = typename Tree::circle;

  int n = wp.size();

  double aveCount = time_loop(
      rounds, 1.0, [&]() {},
      [&]() {
        parlay::parallel_for(0, queryNum, [&](size_t i) {
          Box query_box_seq = pkd.get_box(
              Box(wp[i], wp[i]), Box(wp[(i + n / 2) % n], wp[(i + n / 2) % n]));
          auto d =
              Tree::p2p_distance(wp[i], wp[(i + n / 2) % n], wp[i].get_dim());
          d = static_cast<Coord>(std::sqrt(d));
          circle cl = circle(wp[i], d);
          kdknn[i] = pkd.range_count(cl);
        });
      },
      [&]() {});

  LOG << aveCount << " " << std::flush;

  return;
}

// NOTE: run range query and check the correct
template <typename Point, typename Tree>
void rangeQuery(parlay::sequence<Point> const& wp, Tree& pkd, Typename* kdknn,
                int const& rounds, int const rec_num, int const rec_type,
                int const DIM, parlay::sequence<Point>& Out) {
  // using Points = typename Tree::Points;

  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, true>(rec_num, rec_type, wp, DIM);
  Out.resize(rec_num * max_size);
  size_t step = Out.size() / rec_num;
  // using ref_t = std::reference_wrapper<Point>;
  // parlay::sequence<ref_t> out_ref(Out.size(), std::ref(Out[0]));
  // parlay::sequence<double> preTime( rec_num, 0 );

  double aveQuery = time_loop(
      rounds, 1.0, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeQuery(
              query_box_seq[i].first, Out.cut(i * step, (i + 1) * step));
          kdknn[i] = size;
        });
      },
      [&]() {});

  LOG << "check range query: " << rec_num << " " << rec_type << " " << max_size
      << ENDL;
  for (int i = 0; i < rec_num; i++) {
    assert(std::cmp_equal(kdknn[i], query_box_seq[i].second.size()));
    // LOG << kdknn[i] << " " << query_box_seq[i].second.size() << " "
    //     << query_box_seq[i].first.first << query_box_seq[i].first.second
    //     << ENDL;
    parlay::sort_inplace(
        Out.cut(i * step, i * step + query_box_seq[i].second.size()));
    parlay::sort_inplace(query_box_seq[i].second);
    for (size_t j = 0; j < query_box_seq[i].second.size(); j++) {
      assert(Out[i * step + j] == query_box_seq[i].second.at(j));
      // if (Out[i * step + j] != query_box_seq[i].second.at(j)) LOG << "wrong
      // "; LOG << Out[j] << " " << query_box_seq[i].second.at(j) << ENDL;
    }
  }

  LOG << aveQuery << " " << std::flush;
  return;
}

//* test range count for fix rectangle
template <typename Point, typename Tree>
void rangeCountFix(parlay::sequence<Point> const& WP, Tree& pkd,
                   Typename* kdknn, int const& rounds, int rec_type,
                   int rec_num, int DIM) {
  // using Tree = Tree;
  // using Points = typename Tree::Points;
  // using Box = typename Tree::Box;

  // int n = WP.size();

  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, false>(rec_num, rec_type, WP, DIM);
  parlay::sequence<size_t> vis_nodes(rec_num), gen_box(rec_num),
      full_box(rec_num), skip_box(rec_num);

  double aveCount = time_loop(
      rounds, 1.0, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeCount(query_box_seq[i].first);

          kdknn[i] = size;
          vis_nodes[i] = logger.vis_node_num;
          gen_box[i] = logger.generate_box_num;
          full_box[i] = logger.full_box_num;
          skip_box[i] = logger.skip_box_num;
        });
        // for (int i = 0; i < rec_num; i++) {
        //     kdknn[i] = pkd.RangeCount(query_box_seq[i].first);
        // }
      },
      [&]() {});

  LOG << aveCount << " " << std::flush;
  LOG << parlay::reduce(vis_nodes.cut(0, rec_num)) / rec_num << " "
      << std::flush;
  LOG << parlay::reduce(gen_box.cut(0, rec_num)) / rec_num << " " << std::flush;
  LOG << parlay::reduce(full_box.cut(0, rec_num)) / rec_num << " "
      << std::flush;
  LOG << parlay::reduce(skip_box.cut(0, rec_num)) / rec_num << " "
      << std::flush;

  return;
}
//
// template<typename Point>
// void rangeCountFixWithLog(const parlay::sequence<Point>& WP, BaseTree<Point>&
// pkd, Typename* kdknn, const int& rounds,
//                           int rec_type, int rec_num, int DIM) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     using Box = typename Tree::Box;
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
//         if (query_box_seq[i].second != kdknn[i]) LOG << "wrong" << ENDL;
//         LOG << query_box_seq[i].second << " " << std::scientific << aveQuery
//         << ENDL;
//     }
//
//     return;
// }
//
//* test range query for fix rectangle
template <typename Point, typename Tree>
void rangeQueryFix(parlay::sequence<Point> const& WP, Tree& pkd,
                   Typename* kdknn, int const& rounds,
                   parlay::sequence<Point>& Out, int rec_type, int rec_num,
                   int DIM) {
  // using Tree = Tree;
  // using Points = typename Tree::Points;
  // using Box = typename Tree::Box;

  auto [query_box_seq, max_size] =
      gen_rectangles<Point, Tree, false>(rec_num, rec_type, WP, DIM);
  parlay::sequence<size_t> vis_nodes(rec_num), gen_box(rec_num),
      full_box(rec_num), skip_box(rec_num);
  Out.resize(rec_num * max_size);

  // int n = WP.size();
  size_t step = Out.size() / rec_num;
  // using ref_t = std::reference_wrapper<Point>;
  // parlay::sequence<ref_t> out_ref( Out.size(), std::ref( Out[0] ) );

  double aveQuery = time_loop(
      rounds, 1.0, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = pkd.RangeQuery(
              query_box_seq[i].first, Out.cut(i * step, (i + 1) * step));

          kdknn[i] = size;
          vis_nodes[i] = logger.vis_node_num;
          gen_box[i] = logger.generate_box_num;
          full_box[i] = logger.full_box_num;
          skip_box[i] = logger.skip_box_num;
        });
      },
      [&]() {});

  LOG << aveQuery << " " << std::flush;
  LOG << parlay::reduce(vis_nodes.cut(0, rec_num)) / rec_num << " "
      << std::flush;
  LOG << parlay::reduce(gen_box.cut(0, rec_num)) / rec_num << " " << std::flush;
  LOG << parlay::reduce(full_box.cut(0, rec_num)) / rec_num << " "
      << std::flush;
  LOG << parlay::reduce(skip_box.cut(0, rec_num)) / rec_num << " "
      << std::flush;
  return;
}
//
// template<typename Point>
// void rangeQuerySerialWithLog(const parlay::sequence<Point>& WP,
// BaseTree<Point>& pkd, Typename* kdknn,
//                              const int& rounds, parlay::sequence<Point>& Out,
//                              int rec_type, int rec_num, int DIM) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     using Box = typename Tree::Box;
//
//     auto [query_box_seq, max_size] = gen_rectangles(rec_num, rec_type, WP,
//     DIM); Out.resize(rec_num * max_size);
//
//     size_t step = Out.size() / rec_num;
//     // using ref_t = std::reference_wrapper<Point>;
//     // parlay::sequence<ref_t> out_ref( Out.size(), std::ref( Out[0] ) );
//
//     for (int i = 0; i < rec_num; i++) {
//         double aveQuery = time_loop(
//             rounds, -1.0, [&]() {},
//             [&]() { kdknn[i] = pkd.range_query_serial(query_box_seq[i].first,
//             Out.cut(i * step, (i + 1) * step)); },
//             [&]() {});
//         if (query_box_seq[i].second != kdknn[i]) LOG << "wrong" << ENDL;
//         LOG << query_box_seq[i].second << " " << std::scientific << aveQuery
//         << ENDL;
//         // LOG << query_box_seq[i].second << " " << std::setprecision(7) <<
//         aveQuery << ENDL;
//     }
//
//     return;
// }
//
// template<typename Point>
// void generate_knn(const uint_fast8_t& Dim, const parlay::sequence<Point>& WP,
// const int K, const char* outFile) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     using Coord = typename Point::Coord;
//     using nn_pair = std::pair<Point, Coord>;
//     // using nn_pair = std::pair<std::reference_wrapper<Point>, Coord>;
//     using ID_type = uint;
//
//     size_t n = WP.size();
//
//     Tree pkd;
//     Points wp = Points(n);
//     parlay::copy(WP.cut(0, n), wp.cut(0, n));
//
//     pkd.build(parlay::make_slice(wp), Dim);
//
//     parlay::sequence<nn_pair> Out(K * n);
//     parlay::sequence<kBoundedQueue<Point, nn_pair>> bq =
//         parlay::sequence<kBoundedQueue<Point, nn_pair>>::uninitialized(n);
//     parlay::parallel_for(0, n, [&](size_t i) {
//         bq[i].resize(Out.cut(i * K, i * K + K));
//         bq[i].reset();
//     });
//
//     std::cout << "begin query" << std::endl;
//     parlay::copy(WP.cut(0, n), wp.cut(0, n));
//     node* KDParallelRoot = pkd.get_root();
//     auto bx = pkd.get_root_box();
//     parlay::parallel_for(0, n, [&](size_t i) {
//         size_t visNodeNum = 0;
//         pkd.k_nearest(KDParallelRoot, wp[i], Dim, bq[i], bx, visNodeNum);
//     });
//     std::cout << "finish query" << std::endl;
//
//     std::ofstream ofs(outFile);
//     if (!ofs.is_open()) {
//         throw("file not open");
//         abort();
//     }
//     size_t m = n * K;
//     ofs << "WeightedAdjacencyGraph" << '\n';
//     ofs << n << '\n';
//     ofs << m << '\n';
//     parlay::sequence<uint64_t> offset(n + 1);
//     parlay::parallel_for(0, n + 1, [&](size_t i) { offset[i] = i * K; });
//     // parlay::parallel_for( 0, n, [&]( size_t i ) {
//     //   for ( size_t j = 0; j < K; j++ ) {
//     //     if ( Out[i * K + j].first == wp[i] ) {
//     //       printf( "%d, self-loop\n", i );
//     //       exit( 0 );
//     //     }
//     //   }
//     // } );
//
//     parlay::sequence<ID_type> edge(m);
//     parlay::parallel_for(0, m, [&](size_t i) { edge[i] = Out[i].first.id; });
//     parlay::sequence<double> weight(m);
//     parlay::parallel_for(0, m, [&](size_t i) { weight[i] = Out[i].second; });
//     for (size_t i = 0; i < n; i++) {
//         ofs << offset[i] << '\n';
//     }
//     for (size_t i = 0; i < m; i++) {
//         ofs << edge[i];
//         ofs << "\n";
//         // ofs << edge[i] << '\n';
//     }
//     for (size_t i = 0; i < m; i++) {
//         ofs << weight[i] << '\n';
//     }
//     ofs.close();
//     return;
// }
//
// template<typename Point>
// void insertOsmByTime(const int Dim, const
// parlay::sequence<parlay::sequence<Point>>& node_by_time, const int rounds,
//                      BaseTree<Point>& pkd, const int K, Typename* kdknn) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     using Box = typename Tree::Box;
//
//     pkd.delete_tree();
//     int time_period_num = node_by_time.size();
//     parlay::sequence<Points> wp(time_period_num);
//     for (int i = 0; i < time_period_num; i++) {
//         wp[i].resize(node_by_time[i].size());
//     }
//
//     double ave = time_loop(
//         rounds, 1.0,
//         [&]() {
//             for (int i = 0; i < time_period_num; i++) {
//                 parlay::copy(node_by_time[i], wp[i]);
//             }
//         },
//         [&]() {
//             for (int i = 0; i < time_period_num; i++) {
//                 pkd.BatchInsert(parlay::make_slice(wp[i]), Dim);
//             }
//         },
//         [&]() { pkd.delete_tree(); });
//
//     // NOTE: begin revert
//     for (int i = 0; i < time_period_num; i++) {
//         parlay::copy(node_by_time[i], wp[i]);
//     }
//     LOG << ENDL;
//     for (int i = 0; i < time_period_num; i++) {
//         parlay::internal::timer t;
//         t.reset(), t.start();
//
//         pkd.BatchInsert(parlay::make_slice(wp[i]), Dim);
//
//         t.stop();
//         LOG << wp[i].size() << " " << t.total_time() << " ";
//
//         if (time_period_num < 12) {
//             Points tmp(wp[0].begin(), wp[0].begin() + batchQueryOsmSize);
//             queryKNN(Dim, tmp, rounds, pkd, kdknn, K, true);
//         } else if (i != 0 && (i + 1) % 12 == 0) {
//             Points tmp(batchQueryOsmSize);
//             parlay::copy(parlay::make_slice(wp[0]), tmp.cut(0,
//             wp[0].size())); parlay::copy(parlay::make_slice(wp[1].begin(),
//             wp[1].begin() + batchQueryOsmSize - wp[0].size()),
//                          tmp.cut(wp[0].size(), tmp.size()));
//             queryKNN(Dim, tmp, rounds, pkd, kdknn, K, true);
//         }
//
//         LOG << ENDL;
//     }
//
//     size_t max_deep = 0;
//     LOG << ave << " " << pkd.getMaxTreeDepth(pkd.get_root(), max_deep) << " "
//     << pkd.getAveTreeHeight() << " "
//         << std::flush;
//
//     return;
// }
//
// template<typename Point, int print = 1>
// void incrementalBuildAndQuery(const int Dim, const parlay::sequence<Point>&
// WP, const int rounds, BaseTree<Point>& pkd,
//                               double stepRatio, const
//                               parlay::sequence<Point>& query_points) {
//     using Tree = BaseTree<Point>;
//     using Points = typename Tree::Points;
//     using node = typename Tree::node;
//     size_t n = WP.size();
//     size_t step = n * stepRatio;
//     Points wp = Points::uninitialized(n);
//
//     pkd.delete_tree();
//
//     parlay::copy(WP, wp);
//     size_t l = 0, r = 0;
//
//     size_t batchSize = static_cast<size_t>(WP.size() * knnBatchInbaRatio);
//     Typename* kdknn = new Typename[batchSize];
//     const int k[3] = {1, 5, 100};
//
//     LOG << "begin insert: " << batchSize << ENDL;
//     size_t cnt = 0;
//     while (l < n) {
//         parlay::internal::timer t;
//         t.reset(), t.start();
//         r = std::min(l + step, n);
//         pkd.BatchInsert(wp.cut(l, r), Dim);
//         t.stop();
//
//         // NOTE: print info
//
//         if ((cnt < 50 && cnt % 5 == 0) || (cnt < 100 && cnt % 10 == 0) ||
//         (cnt % 100 == 0)) {
//             size_t max_deep = 0;
//             LOG << l << " " << r << " " << t.total_time() << " " <<
//             pkd.getMaxTreeDepth(pkd.get_root(), max_deep) << " "
//                 << pkd.getAveTreeHeight() << " " << std::flush;
//
//             // NOTE: add additional query phase
//             for (int i = 0; i < 3; i++) {
//                 queryKNN<Point, 0, 1>(Dim, query_points, 1, pkd, kdknn, k[i],
//                 true);
//             }
//             LOG << ENDL;
//         }
//         cnt++;
//         l = r;
//     }
//     delete[] kdknn;
//
//     return;
// }

template <typename T>
class counter_iterator {
 private:
  struct accept_any {
    template <typename U>
    accept_any& operator=(U const&) {
      return *this;
    }
  };

 public:
  typedef std::output_iterator_tag iterator_category;

  counter_iterator(T& counter) : counter(counter) {}
  counter_iterator& operator=(counter_iterator const& other) {
    if (this != &other) {  // Check for self-assignment
      counter = other.counter;
    }
    return *this;
  }
  counter_iterator(counter_iterator const& other) : counter(other.counter) {}

  bool operator==(counter_iterator const& rhs) const {
    return counter == rhs.counter;
  }
  bool operator!=(counter_iterator const& rhs) const {
    return counter != rhs.counter;
  }

  accept_any operator*() const {
    ++counter.get();
    return {};
  }

  counter_iterator& operator++() {  // ++a
    return *this;
  }
  counter_iterator operator++(int) {  // a++
    return *this;
  }

 protected:
  std::reference_wrapper<T> counter;
};

//*---------- generate Points within a 0-box_size --------------------
template <typename Point>
void generate_random_points(parlay::sequence<Point>& wp, Coord _box_size,
                            long n, int Dim) {
  Coord box_size = _box_size;

  std::random_device rd;      // a seed source for the random number engine
  std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distrib(1, box_size);

  parlay::random_generator gen(distrib(gen_mt));
  std::uniform_int_distribution<int> dis(0, box_size);

  wp.resize(n);
  // generate n random Points in a cube
  parlay::parallel_for(
      0, n,
      [&](long i) {
        auto r = gen[i];
        for (int j = 0; j < Dim; j++) {
          wp[i].pnt[j] = dis(r);
        }
      },
      1000);
  return;
}

template <typename Point>
std::pair<size_t, int> read_points(char const* iFile,
                                   parlay::sequence<Point>& wp,
                                   [[maybe_unused]] int K,
                                   [[maybe_unused]] bool withID = false) {
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  static Coords samplePoint;
  parlay::sequence<char> S = readStringFromFile(iFile);
  parlay::sequence<char*> W = stringToWords(S);
  size_t N = std::stoul(W[0], nullptr, 10);
  int Dim = atoi(W[1]);
  assert(N >= 0 && Dim >= 1);

  auto pts = W.cut(2, W.size());
  assert(pts.size() % Dim == 0);
  size_t n = pts.size() / Dim;
  auto a = parlay::tabulate(Dim * n, [&](size_t i) -> Coord {
    if constexpr (std::is_integral_v<Coord>)
      return std::stol(pts[i]);
    else if (std::is_floating_point_v<Coord>)
      return std::stod(pts[i]);
  });
  wp.resize(N);
  parlay::parallel_for(0, n, [&](size_t i) {
    for (int j = 0; j < Dim; j++) {
      wp[i].pnt[j] = a[i * Dim + j];
      if constexpr (std::is_same_v<Point,
                                   PointType<Coord, samplePoint.size()>>) {
        ;
      } else {
        wp[i].id = i;
      }
    }
  });
  return std::make_pair(N, Dim);
}

struct wrapper {
  struct QadTree {
    template <class Point>
    struct Desc {
      using TreeType = cpdd::OrthTree<Point, cpdd::RotateDim<Point>, 2>;
    };
  };
  struct OctTree {
    template <class Point>
    struct Desc {
      using TreeType = cpdd::OrthTree<Point, cpdd::RotateDim<Point>, 3>;
    };
  };
  struct KDtree {
    template <class Point>
    struct Desc {
      using TreeType = cpdd::KdTree<Point, cpdd::MaxStretchDim<Point>>;
      // using TreeType = cpdd::KdTree<Point, cpdd::RotateDim<Point>>;
    };
  };
};
