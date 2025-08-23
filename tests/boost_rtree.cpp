// NOTE: R tree references
// https://www.boost.org/doc/libs/1_86_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
#include <boost/geometry.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iterator>

#include "test_framework.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using Typename = Coord;
// Define a 2D point
// Define a box (for range queries)

template <typename T = void>
class counting_output_iterator {
 private:
  size_t* counter_;

 public:
  // Iterator traits
  using iterator_category = std::output_iterator_tag;
  using value_type = void;
  using difference_type = void;
  using pointer = void;
  using reference = void;

  // Constructor
  explicit counting_output_iterator(size_t& counter) : counter_(&counter) {}

  // Proxy class for assignment
  class proxy {
   private:
    size_t* counter_;

   public:
    explicit proxy(size_t* counter) : counter_(counter) {}

    template <typename U>
    proxy& operator=(U const&) {
      ++(*counter_);
      return *this;
    }
  };

  // Dereference operator returns proxy
  proxy operator*() { return proxy(counter_); }

  // Pre-increment (no-op, but required for output iterator)
  counting_output_iterator& operator++() { return *this; }

  // Post-increment (no-op, but required for output iterator)
  counting_output_iterator operator++(int) { return *this; }
};

// Helper function to create the iterator
template <typename T = void>
counting_output_iterator<T> make_counting_output_iterator(size_t& counter) {
  return counting_output_iterator<T>(counter);
}

// set the points for the R-tree
template <size_t Dim, size_t MaxDim>
struct SetPoints {
  static void set(auto& points, auto& points_insert, auto const& wp,
                  auto const& wi, size_t i) {
    bg::set<Dim>(points[i], wp[i].pnt[Dim]);
    bg::set<Dim>(points_insert[i], wi[i].pnt[Dim]);
    SetPoints<Dim + 1, MaxDim>::set(points, points_insert, wp, wi, i);
  }
};

template <size_t MaxDim>
struct SetPoints<MaxDim, MaxDim> {
  static void set(auto&, auto&, auto const&, auto const&, size_t) {
    // End of recursion
  }
};
template <class TreeDesc, typename RPoint, typename Point, size_t kRTreeMaxEle>
void TestRtreeParallel(int Dim, parlay::sequence<Point>& wp,
                       parlay::sequence<Point>& wi, [[maybe_unused]] int N,
                       int K, [[maybe_unused]] int const& kRounds,
                       int const& kTag, int const& kQueryType,
                       int const kSummary) {}

int main(int argc, char* argv[]) {
  commandLine params(
      argc, argv,
      "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
      "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
      "<_insertFile>] [-s <kSummary>]");

  char* input_file_path = params.getOptionValue("-p");
  int K = params.getOptionIntValue("-k", 10);
  int dims = params.getOptionIntValue("-d", 3);
  size_t N = params.getOptionLongValue("-n", -1);
  int tag = params.getOptionIntValue("-t", 1);
  int rounds = params.getOptionIntValue("-r", 3);
  int query_type = params.getOptionIntValue("-q", 0);
  int read_insert_file = params.getOptionIntValue("-i", 1);
  int summary = params.getOptionIntValue("-s", 0);
  int tree_type = params.getOptionIntValue("-T", 0);
  int split_type = params.getOptionIntValue("-l", 0);

  auto run = [&]<typename TreeDesc, typename Point>(
                 int const& kDim, parlay::sequence<Point> const& wp,
                 parlay::sequence<Point> const& wi, size_t const& N,
                 int const& K, int const& kRounds, string const& kInsertFile,
                 int const& kTag, int const& kQueryType, int const kSummary) {
    // constexpr std::array<size_t, 5> const kMaxEleArr = {2, 4, 8, 32, 100};
    // auto callTestRtreeParallel = [&]<size_t...
    // Is>(std::index_sequence<Is...>) {
    //   (TestRtreeParallel<TreeWrapper,
    //                      bg::model::point<Coord, kDim, bg::cs::cartesian>,
    //                      Point, kMaxEleArr[Is]>(kDim, wp, wi, N, K, rounds,
    //                                             tags, query_type, summary),
    //    ...);
    // };
    // callTestRtreeParallel(std::make_index_sequence<kMaxEleArr.size()>{});
    using Tree = TreeDesc::TreeType;
    using Points = typename Tree::Points;
    using RPoint = bg::model::point<Coord, Point::GetDim(), bg::cs::cartesian>;
    using RBox = bg::model::box<RPoint>;
    constexpr size_t kRTreeMaxEle = 32;  // Max elements per node in R-tree
    using BoostRTree = bgi::rtree<RPoint, bgi::quadratic<kRTreeMaxEle>>;
    // using Points = typename Tree::Points;
    // using Box = typename Tree::Box;

    // Sample points to insert into the R-tree
    // std::cout << kRTreeMaxEle << " ";

    std::vector<RPoint> _points(wp.size());
    // assert(wp.size() == wi.size());

    // NOTE: set the value of points
    auto set_points = [&]<size_t... Is>(auto& r_pt, auto& pt,
                                        std::index_sequence<Is...>) {
      (bg::set<Is>(r_pt, pt.pnt[Is]), ...);
    };

    parlay::parallel_for(0, wp.size(), [&](size_t i) {
      set_points(_points[i], wp[i],
                 std::make_index_sequence<Point::GetDim()>{});
    });

    // To generate box
    auto generate_query_box = [&](int rec_num, int rec_total_type,
                                  Points const& wp) {
      // NOTE: generate rectangles for the first half of the points
      parlay::sequence<parlay::sequence<std::pair<typename Tree::Box, size_t>>>
          query_box_seq(rec_total_type);
      parlay::sequence<size_t> query_max_size(rec_total_type);
      for (int i = 0; i < rec_total_type; i++) {
        auto [query_box, max_size] =
            gen_rectangles<Point, Tree, false, true>(rec_num, i, wp, kDim);
        query_box_seq[i] = query_box;
        query_max_size[i] = max_size;
      }
      return std::make_pair(query_box_seq, query_max_size);
    };

    // summary table
    if (tag & (1 << 3) || tag & (1 << 4)) {
      auto [query_box_seq, query_max_size] =
          generate_query_box(kRangeQueryNum, 3, wp.subseq(0, wp.size() / 2));

      auto run_rtree_knn = [&](BoostRTree& tree, Points const& query_points,
                               int kth) {
        std::vector<std::vector<RPoint>> ans(query_points.size(),
                                             std::vector<RPoint>());
        double ave_knn = time_loop(
            rounds, 1.0, []() {},
            [&]() {
              parlay::parallel_for(0, query_points.size(), [&](size_t i) {
                RPoint query_point(wp[i].pnt[0], wp[i].pnt[1]);
                tree.query(bgi::nearest(query_point, kth),
                           std::back_inserter(ans[i]));
              });
            },
            []() {});
        std::cout << ave_knn << " " << std::flush;
      };

      auto run_rtree_range_count = [&](BoostRTree& tree,
                                       auto const& query_box_seq,
                                       size_t const max_size, int type) {
        size_t query_num = query_box_seq.size();

        double ave_count = time_loop(
            kRounds, -1.0, [&]() {},
            [&]() {
              parlay::parallel_for(0, query_num, [&](size_t s) {
                RPoint a, b;
                set_points(a, query_box_seq[s].first.first,
                           std::make_index_sequence<Point::GetDim()>{});
                set_points(b, query_box_seq[s].first.second,
                           std::make_index_sequence<Point::GetDim()>{});

                RBox query_box(a, b);
                size_t ans = 0;
                tree.query(bgi::within(query_box),
                           make_counting_output_iterator<size_t>(ans));
              });
            },
            [&]() {});
        std::cout << ave_count << " " << std::flush;
      };

      auto run_rtree_range_query = [&](BoostRTree& tree,
                                       auto const& query_box_seq,
                                       size_t const max_size, int type) {
        size_t query_num = query_box_seq.size();
        std::vector<std::vector<RPoint>> ans(query_num, std::vector<RPoint>());

        double aveQuery = time_loop(
            kRounds, -1.0, [&]() {},
            [&]() {
              parlay::parallel_for(0, query_num, [&](size_t s) {
                // RBox query_box(RPoint(queryBox[s].first.first.pnt[0],
                //                       queryBox[s].first.first.pnt[1]),
                //                RPoint(queryBox[s].first.second.pnt[0],
                //                       queryBox[s].first.second.pnt[1]));
                RPoint a, b;
                set_points(a, query_box_seq[s].first.first,
                           std::make_index_sequence<Point::GetDim()>{});
                set_points(b, query_box_seq[s].first.second,
                           std::make_index_sequence<Point::GetDim()>{});

                RBox query_box(a, b);
                tree.query(bgi::within(query_box), std::back_inserter(ans[s]));
              });
            },
            [&]() {});
        std::cout << aveQuery << " " << std::flush;
      };

      auto run_all_tests = [&](BoostRTree& tree) {
        {
          int k[3] = {1, 10, 100};

          std::cout << "in-dis-skewed knn time: ";
          size_t batch_size = static_cast<size_t>(wp.size() * kBatchQueryRatio);
          for (int i = 0; i < 3; i++) {
            run_rtree_knn(tree, wp.subseq(0, batch_size), k[i]);
          }
          puts("");

          std ::cout << "out-dis-skewed knn time: ";
          for (int i = 0; i < 3; i++) {
            run_rtree_knn(tree, wp.subseq(wp.size() - batch_size, wp.size()),
                          k[i]);
          }
          puts("");

          // NOTE: sample points within the whole input datasets
          auto query_pts = parlay::pack(
              wp, parlay::tabulate(wp.size(), [&](size_t i) -> bool {
                return i % (wp.size() / (batch_size * 2)) == 0;
              }));

          std::cout << "in-dis-uniform knn time: ";
          for (int i = 0; i < 3; i++) {
            run_rtree_knn(
                tree, parlay::random_shuffle(query_pts.subseq(0, batch_size)),
                k[i]);
          }
          puts("");

          std::cout << "out-dis-uniform knn time: ";
          for (int i = 0; i < 3; i++) {
            run_rtree_knn(tree,
                          parlay::random_shuffle(
                              query_pts.subseq(batch_size, query_pts.size())),
                          k[i]);
          }
          puts("");
        }

        // NOTE: range count
        {
          std::cout << "range count time: ";
          for (int i = 0; i < 3; i++) {
            run_rtree_range_count(tree, query_box_seq[i], query_max_size[i], i);
          }
          puts("");
        }

        // NOTE: range query
        {
          std::cout << "range query time: ";
          for (int i = 0; i < 3; i++) {
            run_rtree_range_query(tree, query_box_seq[i], query_max_size[i], i);
          }
          puts("");
        }
      };

      // Main test
      std::cout << std::fixed << std::setprecision(5);
      puts("");
      parlay::internal::timer timer;
      BoostRTree tree;

      // incre insert full
      tree = BoostRTree();  // reset tree
      timer.reset(), timer.start();
      for (int i = 0; i < _points.size(); i++) {
        tree.insert(_points[i]);
      }
      std::cout << "## incre insert full: " << timer.total_time() << " "
                << std::endl;
      run_all_tests(tree);

      // incre insert half
      tree = BoostRTree();  // reset tree
      timer.reset(), timer.start();
      for (int i = 0; i < _points.size() / 2; i++) {
        tree.insert(_points[i]);
      }
      std::cout << "## incre insert half: " << timer.total_time() << " "
                << std::endl;
      run_all_tests(tree);

      // directly build half
      timer.reset(), timer.start();
      tree = BoostRTree(_points.begin(), _points.begin() + _points.size() / 2);
      std::cout << "## build tree half: " << timer.total_time() << " "
                << std::endl;
      run_all_tests(tree);

      tree = BoostRTree(_points.begin(), _points.end());
      timer.reset(), timer.start();
      for (int i = _points.size() / 2; i < _points.size(); i++) {
        tree.remove(_points[i]);
      }
      std::cout << "## incre delete half: " << timer.total_time() << " "
                << std::endl;
      run_all_tests(tree);
    }

    if (tag & (1 << 6)) {
      puts("");
      // incre insert half
      parlay::internal::timer timer;
      auto tree = BoostRTree();  // reset tree
      timer.reset(), timer.start();
      for (int i = 0; i < _points.size() / 2; i++) {
        tree.insert(_points[i]);
      }
      std::cout << "## incre insert half: " << timer.total_time() << " "
                << std::endl;

      auto [query_box_seq, query_max_size] = generate_query_box(
          kSingleQueryLogRepeatNum, 3, wp.subseq(0, wp.size() / 2));
      int rec_num = kSingleQueryLogRepeatNum;
      size_t query_num = rec_num;
      std::vector<std::vector<RPoint>> ans(query_num, std::vector<RPoint>());

      auto run_range_query_log = [&](int rec_type, auto const& query_box_seq,
                                     auto query_max_size) {
        for (int i = 0; i < rec_num; i++) {
          parlay::internal::timer t;
          t.reset(), t.start();
          RPoint a, b;
          set_points(a, query_box_seq[i].first.first,
                     std::make_index_sequence<Point::GetDim()>{});
          set_points(b, query_box_seq[i].first.second,
                     std::make_index_sequence<Point::GetDim()>{});

          RBox query_box(a, b);
          tree.query(bgi::within(query_box), std::back_inserter(ans[i]));
          t.stop();
          std::cout << rec_type << " " << query_box_seq[i].second << " "
                    << std::scientific << t.total_time() << std::endl;
        }
      };

      for (int rec_type = 0; rec_type < 3; rec_type++) {
        run_range_query_log(rec_type, query_box_seq[rec_type],
                            query_max_size[rec_type]);
      }
    }

    // if (kTag & (1 << 0)) {  // insert
    //   auto rtree_insert = [&](auto r_tree, double r) {
    //     timer.reset();
    //     timer.start();
    //     size_t sz = _points_insert.size() * r;
    //     r_tree.insert(_points_insert.begin(), _points_insert.begin() + sz);
    //     std::cout << timer.total_time() << " " << std::flush;
    //   };

    //   if (kSummary) {
    //     parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
    //     for (size_t i = 0; i < ratios.size(); i++) {
    //       auto r_tree = tree;
    //       rtree_insert(r_tree, ratios[i]);
    //     }
    //   } else {
    //     auto r_tree = tree;
    //     rtree_insert(r_tree, kBatchInsertRatio);
    //   }

    //   // if (kTag == 1) wp.append(wi);
    // }

    // if (kTag & (1 << 1)) {  // delete
    //   auto rtree_delete = [&](auto& r_tree, double ratio = 1.0) {
    //     timer.reset();
    //     timer.start();
    //     assert(tree.size() == wp.size());
    //     size_t sz = _points.size() * ratio;
    //     r_tree.remove(_points.begin(), _points.begin() + sz);
    //     timer.stop();
    //     std::cout << timer.total_time() << " " << std::flush;
    //   };

    //   if (kSummary) {
    //     parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
    //     for (size_t i = 0; i < ratios.size(); i++) {
    //       auto r_tree = tree;
    //       rtree_delete(r_tree, ratios[i]);
    //     }
    //   } else {
    //     auto r_tree = tree;
    //     rtree_delete(r_tree, kBatchInsertRatio);
    //   }
    // }

    // kNN query: find the 3 nearest neighbors to the point (2.5, 2.5)
    // Point query_point(2.5, 2.5);
    // std::vector<Point> knn_results;
    // rtree.query(bgi::nearest(query_point, 3),
    // std::back_inserter(knn_results));
    //
    // std::cout << "kNN query results (3 nearest to (2.5, 2.5)):" <<
    // std::std::endl; for (const auto& v : knn_results) {
    //     std::cout << "Point: (" << v.first.get<0>() << ", " <<
    //     v.first.get<1>()
    //     << "), ID: " << v.second << std::std::endl;
    // }
    // if (kQueryType & (1 << 0)) {  // NOTE: KNN query
    //   auto run_rtree_knn = [&](int kth, size_t batchSize) {
    //     timer.reset();
    //     timer.start();
    //     parlay::sequence<size_t> visNodeNum(batchSize, 0);
    //     parlay::parallel_for(0, batchSize, [&](size_t i) {
    //       RPoint query_point(wp[i].pnt[0], wp[i].pnt[1]);
    //       std::vector<RPoint> knn_results;
    //       tree.query(bgi::nearest(query_point, kth),
    //                  std::back_inserter(knn_results));
    //     });
    //     timer.stop();
    //     std::cout << timer.total_time() << " " << std::flush;
    //   };
    //
    //   size_t batchSize = static_cast<size_t>(wp.size() * kBatchQueryRatio);
    //   if (kSummary == 0) {
    //     int const k[3] = {1, 10, 100};
    //     for (int i = 0; i < 3; i++) {
    //       run_rtree_knn(k[i], batchSize);
    //     }
    //   } else {
    //     run_rtree_knn(K, batchSize);
    //   }
    // }

    //
    // // Range query: find points within the box defined by (1.5, 1.5) and
    // (4.5, 4.5) Box query_box(Point(1.5, 1.5), Point(4.5, 4.5));
    // std::vector<Point> range_results;
    // rtree.query(bgi::intersects(query_box),
    // std::back_inserter(range_results));
    //
    // if (kQueryType & (1 << 2)) {  // NOTE: range query
    //   auto run_rtree_range_query = [&](int type) {
    //     int queryNum = kSummary ? kSummaryRangeQueryNum : kRangeQueryNum;
    //     auto [queryBox, maxSize] =
    //         gen_rectangles<Point, Tree, false, true>(queryNum, type, wp,
    //         kDim);
    //     // using ref_t = std::reference_wrapper<Point_d>;
    //     // std::vector<ref_t> out_ref( queryNum * maxSize, std::ref(
    //     _points[0]
    //     // )
    //     // );
    //     std::vector<RPoint> _ans(queryNum * maxSize);
    //
    //     double aveQuery = time_loop(
    //         kRounds, -1.0, [&]() {},
    //         [&]() {
    //           parlay::parallel_for(0, queryNum, [&](size_t s) {
    //             // RBox query_box(RPoint(queryBox[s].first.first.pnt[0],
    //             //                       queryBox[s].first.first.pnt[1]),
    //             //                RPoint(queryBox[s].first.second.pnt[0],
    //             //                       queryBox[s].first.second.pnt[1]));
    //             RPoint a, b;
    //             set_points(a, queryBox[s].first.first,
    //                        std::make_index_sequence<Point::GetDim()>{});
    //             set_points(b, queryBox[s].first.second,
    //                        std::make_index_sequence<Point::GetDim()>{});
    //
    //             RBox query_box(a, b);
    //             std::vector<RPoint> range_results;
    //             tree.query(bgi::within(query_box),
    //                        std::back_inserter(range_results));
    //           });
    //         },
    //         [&]() {});
    //     std::cout << aveQuery << " " << std::flush;
    //   };
    //
    //   if (kSummary == 0) {
    //     int const type[3] = {0, 1, 2};
    //     for (int i = 0; i < 3; i++) {
    //       run_rtree_range_query(type[i]);
    //     }
    //   } else {
    //     run_rtree_range_query(2);
    //   }
    // }
    //
    // // Range query: find points within the box defined by (1.5, 1.5) and
    // (4.5, 4.5) Box query_box(Point(1.5, 1.5), Point(4.5, 4.5));
    // std::vector<Point> range_results;
    // rtree.query(bgi::intersects(query_box),
    // std::back_inserter(range_results));
    //
    // std::cout << "\nRange query results (points within box (1.5, 1.5) to
    // (4.5, 4.5)):" << std::std::endl; for (const auto& v : range_results) {
    //     std::cout << "Point: (" << v.first.get<0>() << ", " <<
    //     v.first.get<1>()
    //     << "), ID: " << v.second << std::std::endl;
    // }
    //
    // // Batch deletion: delete points within the box (1.5, 1.5) to (4.5, 4.5)
    // for (const auto& v : range_results) {
    //     rtree.remove(v);
    // }
    //
    // // Confirm deletion with another range query
    // std::vector<Point> post_delete_results;
    // rtree.query(bgi::intersects(query_box),
    // std::back_inserter(post_delete_results));
    //
    // std::cout << "\nRange query results after deletion:" << std::std::endl;
    // if (post_delete_results.empty()) {
    //     std::cout << "No points found within the box (1.5, 1.5) to
    //     (4.5, 4.5)"
    //     << std::std::endl;
    // } else {
    //     for (const auto& v : post_delete_results) {
    //         std::cout << "Point: (" << v.first.get<0>() << ", " <<
    //         v.first.get<1>() << "), ID: " << v.second
    //                   << std::std::endl;
    //     }
    // }
    std::cout << std::endl;
    return;
  };

  Wrapper::ApplyOrthogonal(tree_type, dims, split_type, params, run);
  return 0;
}
