// NOTE: R tree references
// https://www.boost.org/doc/libs/1_86_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <vector>

#include "common/parse_command_line.h"
#include "parlay/internal/get_time.h"
#include "parlay/parallel.h"
#include "parlay/slice.h"
#include "test_framework.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using Typename = Coord;
// Define a 2D point
typedef bg::model::point<Coord, 2, bg::cs::cartesian> RPoint;
// Define a box (for range queries)
typedef bg::model::box<RPoint> RBox;

template <class TreeDesc, typename Point, size_t kRTreeMaxEle>
void TestRtreeParallel(int Dim, parlay::sequence<Point>& wp,
                       parlay::sequence<Point>& wi, [[maybe_unused]] int N,
                       int K, [[maybe_unused]] int const& kRounds,
                       int const& kTag, int const& kQueryType,
                       int const kSummary) {
  using Tree = TreeDesc::TreeType;
  // using Points = typename Tree::Points;
  // using Box = typename Tree::Box;

  // Sample points to insert into the R-tree
  std::cout << kRTreeMaxEle << " ";

  std::vector<RPoint> _points(wp.size());
  std::vector<RPoint> _points_insert(wi.size());
  // assert(wp.size() == wi.size());
  parlay::parallel_for(0, wi.size(), [&](size_t i) {
    bg::set<0>(_points[i], wp[i].pnt[0]);
    bg::set<0>(_points_insert[i], wi[i].pnt[0]);

    bg::set<1>(_points[i], wp[i].pnt[1]);
    bg::set<1>(_points_insert[i], wi[i].pnt[1]);
  });

  parlay::internal::timer timer;
  timer.start();
  bgi::rtree<RPoint, bgi::rstar<kRTreeMaxEle>> tree(_points.begin(),
                                                    _points.end());
  timer.stop();
  std::cout << timer.total_time() << " " << -1 << " " << std::flush;

  if (kTag & (1 << 0)) {
    auto rtree_insert = [&](double r) {
      timer.reset();
      timer.start();
      size_t sz = _points_insert.size() * r;
      tree.insert(_points_insert.begin(), _points_insert.begin() + sz);
      std::cout << timer.total_time() << " " << std::flush;
    };

    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        tree.clear();
        tree.insert(_points.begin(), _points.end());
        rtree_insert(ratios[i]);
      }
    } else {
      rtree_insert(kBatchInsertRatio);
    }

    if (kTag == 1) wp.append(wi);
  }

  if (kTag & (1 << 1)) {
    auto rtree_delete = [&](bool afterInsert = 1, double ratio = 1.0) {
      if (!afterInsert) {
        tree.clear();
        tree.insert(_points.begin(), _points.end());
      }
      timer.reset();
      timer.start();
      if (afterInsert) {
        size_t sz = _points_insert.size() * ratio;
        tree.remove(_points_insert.begin(), _points_insert.begin() + sz);
      } else {
        assert(tree.size() == wp.size());
        size_t sz = _points.size() * ratio;
        tree.remove(_points.begin(), _points.begin() + sz);
      }
      timer.stop();
      std::cout << timer.total_time() << " " << std::flush;
    };

    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        rtree_delete(0, ratios[i]);
      }
      tree.clear();
      tree.insert(_points.begin(), _points.end());
    } else {
      rtree_delete(0, kBatchInsertRatio);
    }
  }

  // PERF: handle the size of cgknn dynamically

  // kNN query: find the 3 nearest neighbors to the point (2.5, 2.5)
  // Point query_point(2.5, 2.5);
  // std::vector<Point> knn_results;
  // rtree.query(bgi::nearest(query_point, 3), std::back_inserter(knn_results));
  //
  // std::cout << "kNN query results (3 nearest to (2.5, 2.5)):" <<
  // std::std::endl; for (const auto& v : knn_results) {
  //     std::cout << "Point: (" << v.first.get<0>() << ", " << v.first.get<1>()
  //     << "), ID: " << v.second << std::std::endl;
  // }
  if (kQueryType & (1 << 0)) {  // NOTE: KNN query
    auto run_rtree_knn = [&](int kth, size_t batchSize) {
      timer.reset();
      timer.start();
      parlay::sequence<size_t> visNodeNum(batchSize, 0);
      parlay::parallel_for(0, batchSize, [&](size_t i) {
        RPoint query_point(wp[i].pnt[0], wp[i].pnt[1]);
        std::vector<RPoint> knn_results;
        tree.query(bgi::nearest(query_point, kth),
                   std::back_inserter(knn_results));
      });
      timer.stop();
      std::cout << timer.total_time() << " " << std::flush;
    };

    size_t batchSize = static_cast<size_t>(wp.size() * batchQueryRatio);
    if (kSummary == 0) {
      int const k[3] = {1, 10, 100};
      for (int i = 0; i < 3; i++) {
        run_rtree_knn(k[i], batchSize);
      }
    } else {
      run_rtree_knn(K, batchSize);
    }
  }

  //
  // // Range query: find points within the box defined by (1.5, 1.5) and
  // (4.5, 4.5) Box query_box(Point(1.5, 1.5), Point(4.5, 4.5));
  // std::vector<Point> range_results;
  // rtree.query(bgi::intersects(query_box), std::back_inserter(range_results));
  //
  if (kQueryType & (1 << 2)) {  // NOTE: range query
    auto run_rtree_range_query = [&](int type) {
      int queryNum = kSummary ? kSummaryRangeQueryNum : kRangeQueryNum;
      auto [queryBox, maxSize] =
          gen_rectangles<Point, Tree, false>(queryNum, type, wp, Dim);
      // using ref_t = std::reference_wrapper<Point_d>;
      // std::vector<ref_t> out_ref( queryNum * maxSize, std::ref( _points[0] )
      // );
      std::vector<RPoint> _ans(queryNum * maxSize);

      double aveQuery = time_loop(
          kRounds, -1.0, [&]() {},
          [&]() {
            parlay::parallel_for(0, queryNum, [&](size_t s) {
              RBox query_box(RPoint(queryBox[s].first.first.pnt[0],
                                    queryBox[s].first.first.pnt[1]),
                             RPoint(queryBox[s].first.second.pnt[0],
                                    queryBox[s].first.second.pnt[1]));
              std::vector<RPoint> range_results;
              tree.query(bgi::within(query_box),
                         std::back_inserter(range_results));
            });
          },
          [&]() {});
      std::cout << aveQuery << " " << std::flush;
    };

    if (kSummary == 0) {
      int const type[3] = {0, 1, 2};
      for (int i = 0; i < 3; i++) {
        run_rtree_range_query(type[i]);
      }
    } else {
      run_rtree_range_query(2);
    }
  }
  //
  // // Range query: find points within the box defined by (1.5, 1.5) and
  // (4.5, 4.5) Box query_box(Point(1.5, 1.5), Point(4.5, 4.5));
  // std::vector<Point> range_results;
  // rtree.query(bgi::intersects(query_box), std::back_inserter(range_results));
  //
  // std::cout << "\nRange query results (points within box (1.5, 1.5) to
  // (4.5, 4.5)):" << std::std::endl; for (const auto& v : range_results) {
  //     std::cout << "Point: (" << v.first.get<0>() << ", " << v.first.get<1>()
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
  //     std::cout << "No points found within the box (1.5, 1.5) to (4.5, 4.5)"
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
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>] [-s <kSummary>]");

  char* input_file_path = P.getOptionValue("-p");
  int K = P.getOptionIntValue("-k", 100);
  int kDim = P.getOptionIntValue("-d", 3);
  size_t N = P.getOptionLongValue("-n", -1);
  int kTag = P.getOptionIntValue("-t", 1);
  int kRounds = P.getOptionIntValue("-r", 3);
  int kQueryType = P.getOptionIntValue("-q", 0);
  int read_insert_file = P.getOptionIntValue("-i", 1);
  int kSummary = P.getOptionIntValue("-s", 0);
  // int tree_type = P.getOptionIntValue("-T", 0);

  auto run_test = [&]<class Wrapper>(Wrapper) {
    auto run = [&](auto dim_wrapper) {
      constexpr auto const kDim = decltype(dim_wrapper)::value;
      using PointTypeAlias = PointType<Coord, kDim>;
      using Points = parlay::sequence<PointTypeAlias>;
      using Desc = typename Wrapper::template Desc<PointTypeAlias>;

      std::string name, insert_file_path = "";
      Points wp, wi;

      if (input_file_path != NULL) {  // NOTE: read main Points
        name = std::string(input_file_path);
        name = name.substr(name.rfind("/") + 1);
        std::cout << name << " \n";
        auto [n, d] = read_points<PointTypeAlias>(input_file_path, wp, K);
        N = n;
        assert(d == kDim);
      }

      if (read_insert_file == 1) {  // NOTE: read Points to be inserted
        int id = std::stoi(name.substr(0, name.find_first_of('.')));
        id = (id + 1) % 3;  // WARN: MOD graph number used to test
        if (!id) id++;
        int pos = std::string(input_file_path).rfind("/") + 1;
        insert_file_path = std::string(input_file_path).substr(0, pos) +
                           std::to_string(id) + ".in";
        [[maybe_unused]] auto [n, d] =
            read_points<PointTypeAlias>(insert_file_path.c_str(), wi, K);
        assert(d == kDim);
      }

      constexpr std::array<size_t, 5> const kMaxEleArr = {2, 4, 8, 32, 100};
      auto callTestRtreeParallel =
          [&]<size_t... Is>(std::index_sequence<Is...>) {
            (TestRtreeParallel<Desc, PointTypeAlias, kMaxEleArr[Is]>(
                 kDim, wp, wi, N, K, kRounds, kTag, kQueryType, kSummary),
             ...);
          };
      callTestRtreeParallel(std::make_index_sequence<kMaxEleArr.size()>{});
      // TestRtreeParallel<Desc>(
      //     kDim, wp, wi, N, K, kRounds, kTag, kQueryType, kSummary);
    };

    if (kTag == -1) {
      // NOTE: serial run
      ;
    } else if (kDim == 2) {
      run(std::integral_constant<int, 2>{});
    } else if (kDim == 3) {
      run(std::integral_constant<int, 3>{});
    }
    // } else if (kDim == 5) {
    //     run(std::integral_constant<int, 5>{});
    // } else if (kDim == 7) {
    //     run(std::integral_constant<int, 7>{});
    // } else if (kDim == 9) {
    //     run(std::integral_constant<int, 9>{});
    // } else if (kDim == 10) {
    //     run(std::integral_constant<int, 10>{});
    // }
  };

  run_test(wrapper::KDtree{});  // providing the data types
  return 0;
}
