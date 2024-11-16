
#include <algorithm>
#include <cstddef>
#include <cstdlib>

#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "parlay/primitives.h"
#include "test_framework.h"

template <class TreeDesc, typename Point>
void TestSpacialTree([[maybe_unused]] int const& kDim,
                     parlay::sequence<Point> const& wp,
                     [[maybe_unused]] parlay::sequence<Point> const& wi,
                     [[maybe_unused]] size_t const& N,
                     [[maybe_unused]] int const& K, int const& kRounds,
                     [[maybe_unused]] string const& kInsertFile,
                     [[maybe_unused]] int const& kTag,
                     [[maybe_unused]] int const& kQueryType,
                     [[maybe_unused]] int const kSummary) {
  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;

  Tree tree;
  constexpr bool kTestTime = true;

  BuildTree<Point, Tree, kTestTime, 2>(wp, kRounds, tree);

  Typename* kdknn = nullptr;

  // // NOTE: batch insert
  // if (kTag & (1 << 0)) {
  //   if (kSummary) {
  //     parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
  //     for (size_t i = 0; i < ratios.size(); i++) {
  //       BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds,
  //       ratios[i]);
  //     }
  //   } else {
  //     BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, 0.1);
  //   }
  // }
  //
  // // NOTE: batch delete
  // if (kTag & (1 << 1)) {
  //   if (kSummary) {
  //     parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
  //     for (size_t i = 0; i < ratios.size(); i++) {
  //       BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds,
  //       ratios[i]);
  //     }
  //   } else {
  //     BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds,
  //                                         kBatchInsertRatio);
  //   }
  // }
  //
  // // NOTE: batch diff
  // if (kTag & (1 << 2)) {
  //   if (kSummary) {
  //     parlay::sequence<double> const overlap_ratios = {0.1, 0.2, 0.5, 1};
  //     for (size_t i = 0; i < overlap_ratios.size(); i++) {
  //       BatchDiff<Point, Tree, kTestTime>(
  //           tree, wp, kRounds, kBatchDiffTotalRatio, overlap_ratios[i]);
  //     }
  //   } else {
  //     BatchDiff<Point, Tree, kTestTime>(tree, wp, kRounds,
  //     kBatchDiffTotalRatio,
  //                                       kBatchDiffOverlapRatio);
  //   }
  // }
  //
  // if (kQueryType & (1 << 0)) {  // NOTE: KNN
  //   auto run_batch_knn = [&](Points const& pts, int kth, size_t batchSize) {
  //     Points newPts(batchSize);
  //     parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
  //     kdknn = new Typename[batchSize];
  //     queryKNN<Point>(kDim, newPts, kRounds, tree, kdknn, kth, true);
  //     delete[] kdknn;
  //   };
  //
  //   size_t batchSize = static_cast<size_t>(wp.size() * batchQueryRatio);
  //
  //   if (kSummary == 0) {
  //     int k[3] = {1, 10, 100};
  //     for (int i = 0; i < 3; i++) {
  //       run_batch_knn(wp, k[i], batchSize);
  //     }
  //   } else {  // test kSummary
  //     run_batch_knn(wp, K, batchSize);
  //   }
  // }
  //
  if (kQueryType & (1 << 1)) {  // NOTE: range count
    int recNum = kRangeQueryNum;
    kdknn = new Typename[recNum];

    // std::cout << std::endl;
    for (int i = 0; i < 3; i++) {
      rangeCountFix<Point>(wp, tree, kdknn, kRounds, i, recNum, kDim);
    }

    delete[] kdknn;
  }

  if (kQueryType & (1 << 2)) {  // NOTE: range query
    if (kSummary == 0) {
      int recNum = kRangeQueryNum;
      kdknn = new Typename[recNum];

      for (int i = 0; i < 3; i++) {
        Points Out;
        rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, i, recNum, kDim);
      }
      delete[] kdknn;
    } else if (kSummary == 1) {  // NOTE: for kSummary
      kdknn = new Typename[kSummaryRangeQueryNum];
      Points Out;
      rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, 2,
                           kSummaryRangeQueryNum, kDim);
      delete[] kdknn;
    }
  }

  std::cout << std::endl << std::flush;

  tree.DeleteTree();

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
  int tree_type = P.getOptionIntValue("-T", 0);

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
        std::cout << name << " ";
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

      TestSpacialTree<Desc>(kDim, wp, wi, N, K, kRounds, insert_file_path, kTag,
                            kQueryType, kSummary);
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

  if (tree_type == 0) {
    run_test(wrapper::KDtree{});
  } else if (tree_type == 1 && kDim == 2) {
    run_test(wrapper::QuadTree{});
  } else if (tree_type == 1 && kDim == 3) {
    run_test(wrapper::OctTree{});
  } else if (tree_type == 2) {
    run_test(wrapper::RTree{});
  }

  return 0;
}
