
#include <algorithm>
#include <cstddef>
#include <cstdlib>

#include "parlay/primitives.h"
#include "pspt/dependence/splitter.h"
#include "pspt/kd_tree.h"
#include "pspt/p_tree_impl/cpamtree.hpp"
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
  // // WARN: compress the kdnode to MultiNode, should remove except for exp
  // // if constexpr (IsKdTree<Tree>) {
  // //   tree.Compress2Multi();
  // // }
  //

  Typename* kdknn = nullptr;
  if (kQueryType & (1 << 0)) {  // NOTE: KNN
    auto run_batch_knn = [&](Points const& pts, int kth, size_t batchSize) {
      Points newPts(batchSize);
      parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
      kdknn = new Typename[batchSize];
      queryKNN<Point>(kDim, newPts, kRounds, tree, kdknn, kth, true);
      delete[] kdknn;
    };

    size_t batchSize = static_cast<size_t>(wp.size() * kBatchQueryRatio);

    if (kSummary == 0) {
      int k[3] = {1, 10, 100};
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp, k[i], batchSize);
      }
    } else {  // test kSummary
      run_batch_knn(wp, K, batchSize);
    }
  }
  //
  // if (kQueryType & (1 << 1)) {  // NOTE: range count
  //   if (!kSummary) {
  //     int recNum = kRangeQueryNum;
  //     kdknn = new Typename[recNum];
  //
  //     // std::cout << std::endl;
  //     for (int i = 0; i < 3; i++) {
  //       rangeCountFix<Point>(wp, tree, kdknn, kRounds, i, recNum, kDim);
  //     }
  //
  //     delete[] kdknn;
  //   }
  // }
  //
  // if (kQueryType & (1 << 2)) {  // NOTE: range query
  //   if (kSummary == 0) {
  //     int recNum = kRangeQueryNum;
  //     kdknn = new Typename[recNum];
  //
  //     for (int i = 0; i < 3; i++) {
  //       Points Out;
  //       rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, i, recNum, kDim);
  //     }
  //     delete[] kdknn;
  //   } else if (kSummary == 1) {  // NOTE: for kSummary
  //     kdknn = new Typename[kSummaryRangeQueryNum];
  //     Points Out;
  //     rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, 2,
  //                          kSummaryRangeQueryNum, kDim);
  //     delete[] kdknn;
  //   }
  // }

  std::cout << "\n" << std::flush;

  tree.DeleteTree();

  return;
}

template <class TreeDesc, typename Point>
void TestCPAMBB([[maybe_unused]] int const& kDim,
                parlay::sequence<Point> const& wp,
                [[maybe_unused]] parlay::sequence<Point> const& wi,
                [[maybe_unused]] size_t const& N, [[maybe_unused]] int const& K,
                int const& kRounds, [[maybe_unused]] string const& kInsertFile,
                [[maybe_unused]] int const& kTag,
                [[maybe_unused]] int const& kQueryType,
                [[maybe_unused]] int const kSummary) {
  // auto P = parlay::sequence<geobase::Point>::uninitialized(wp.size());
  // parlay::parallel_for(0, wp.size(), [&](int i){
  //   P[i].x = wp[i].pnt[0];
  //   P[i].y = wp[i].pnt[1];
  //   P[i].id = i;
  // });
  // auto tree = CPAMTree::map_init(P, true);
  // cout << tree.size() << endl;
  auto ptree = BuildPTree(wp, kRounds);
  // test sort time
  return;

  // NOTE: batch insert
  if (kTag & (1 << 0)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchInsertPTree(ptree, wp, wi, kRounds, ratios[i]);
      }
    } else {
      BatchInsertPTree(ptree, wp, wi, kRounds, 0.1);
    }
  }

  // NOTE: batch delete
  if (kTag & (1 << 1)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchDeletePTree(ptree, wp, wi, kRounds, ratios[i]);
      }
    } else {
      BatchDeletePTree(ptree, wp, wi, kRounds, 1.0);
    }
  }

  if (kTag & (1 << 2)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        std::cout << "-1 " << std::flush;
        // BatchDeletePTree(ptree, wp, wi, kRounds, ratios[i]);
      }
    } else {
      std::cout << "-1 " << std::flush;
      // BatchDeletePTree(ptree, wp, wi, kRounds, 1.0);
    }
  }

  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;
  Typename* kdknn = nullptr;

  // Tree kdtree;
  // constexpr bool kTestTime = true;

  // BuildTree<Point, Tree, kTestTime, 2>(wp, kRounds, kdtree);
  if (kQueryType & (1 << 0)) {  // NOTE: KNN
    auto run_batch_knn = [&](Points const& pts, int kth, size_t batchSize) {
      std::cout << "-1 -1 -1 -1 -1 " << std::flush;
      // Points newPts(batchSize);
      // parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
      // kdknn = new Typename[batchSize];
      // queryKNN<Point>(kDim, newPts, kRounds, tree, kdknn, kth, true);
      // delete[] kdknn;
    };

    size_t batchSize = static_cast<size_t>(wp.size() * kBatchQueryRatio);

    if (kSummary == 0) {
      int k[3] = {1, 10, 100};
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp, k[i], batchSize);
      }
    } else {  // test kSummary
      run_batch_knn(wp, K, batchSize);
    }
  }
  
  if (kQueryType & (1 << 1)) {  // NOTE: range count
    if (!kSummary) {
      int recNum = kRangeQueryNum;
      kdknn = new Typename[recNum];

      for (int i = 0; i < 3; i++) {
        rangeCountPtree<Point, Tree>(wp, ptree, kdknn, kRounds, i, recNum,
                                     kDim);
      }

      delete[] kdknn;
    }
  }

  if (kQueryType & (1 << 2)) {  // NOTE: range query
    if (kSummary == 0) {
      int recNum = kRangeQueryNum;
      kdknn = new Typename[recNum];

      for (int i = 0; i < 3; i++) {
        Points Out;
        // rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, i, recNum, kDim);
        rangeReportPtree<Point, Tree>(wp, ptree, kdknn, kRounds, i, recNum,
                                      kDim);
      }
      delete[] kdknn;
    } else if (kSummary == 1) {  // NOTE: for kSummary
      kdknn = new Typename[kSummaryRangeQueryNum];
      Points Out;
      // rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, 2,
      //                      kSummaryRangeQueryNum, kDim);
      rangeReportPtree<Point, Tree>(wp, ptree, kdknn, kRounds, 2,
                                    kSummaryRangeQueryNum, kDim);
      delete[] kdknn;
    }
  }

  // if (kQueryType & (1 << 2)) {  // NOTE: range count
  //   if (!kSummary) {
  //     int recNum = kRangeQueryNum;
  //     kdknn = new Typename[recNum];
  //
  //     // std::cout << std::endl;
  //     for (int i = 0; i < 3; i++) {
  //       rangeReportPtree<Point, Tree>(wp, ptree, kdknn, kRounds, i, recNum,
  //                                     kDim);
  //     }
  //
  //     delete[] kdknn;
  //   }
  // }
  // auto P_insert = P.substr(0, 1234);
  // parlay::parallel_for(0, P_insert.size(), [&](int i){
  //   P_insert[i].id += P.size();
  // });
  // auto tree2 = CPAMTree::map_insert(P_insert, tree, true);
  // cout << tree2.size() << endl;

  // auto P_delete = P_insert.substr(0, 234);
  // auto tree3 = CPAMTree::map_delete(P_delete, tree2, true);
  // cout << tree3.size() << endl;
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>] [-s <kSummary>]");

  char* input_file_path = P.getOptionValue("-p");
  int K = P.getOptionIntValue("-k", 100);
  int dims = P.getOptionIntValue("-d", 3);
  size_t N = P.getOptionLongValue("-n", -1);
  int tag = P.getOptionIntValue("-t", 1);
  int rounds = P.getOptionIntValue("-r", 3);
  int query_type = P.getOptionIntValue("-q", 0);
  int read_insert_file = P.getOptionIntValue("-i", 1);
  int summary = P.getOptionIntValue("-s", 0);
  int tree_type = P.getOptionIntValue("-T", 0);
  int split_type = P.getOptionIntValue("-l", 0);

  auto run = [&]<typename TreeWrapper>() {
    using Point = typename TreeWrapper::Point;
    using Points = parlay::sequence<Point>;
    constexpr auto kDim = Point::GetDim();

    PrintTreeParam<TreeWrapper>();

    std::string name, insert_file_path = "";
    Points wp, wi;

    if (input_file_path != NULL) {  // NOTE: read main Points
      name = std::string(input_file_path);
      name = name.substr(name.rfind('/') + 1);
      std::cout << name << " ";
      auto [n, d] = read_points<Point>(input_file_path, wp, K);
      N = n;
      assert(d == kDim);
    }

    if (read_insert_file == 1) {  // NOTE: read Points to be inserted
      int id = std::stoi(name.substr(0, name.find_first_of('.')));
      id = (id + 1) % 3;  // WARN: MOD graph number used to test
      if (!id) id++;
      auto pos = std::string(input_file_path).rfind('/') + 1;
      insert_file_path = std::string(input_file_path).substr(0, pos) +
                         std::to_string(id) + ".in";
      [[maybe_unused]] auto [n, d] =
          read_points<Point>(insert_file_path.c_str(), wi, K);
      assert(d == kDim);
    }

    TestCPAMBB<TreeWrapper, Point>(kDim, wp, wi, N, K, rounds, insert_file_path,
                                   tag, query_type, summary);

    // TestSpacialTree<TreeWrapper, Point>(
    //     kDim, wp, wi, N, K, rounds, insert_file_path, tag, query_type,
    //     summary);
  };

  Wrapper::Apply(tree_type, dims, split_type, run);

  return 0;
}
