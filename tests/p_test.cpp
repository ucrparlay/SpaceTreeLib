#include "test_framework.h"

int main(int argc, char* argv[]) {
  commandLine params(
      argc, argv,
      "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
      "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
      "<_insertFile>] [-s <kSummary>]");

  char* input_file_path = params.getOptionValue("-p");
  int K = params.getOptionIntValue("-k", 100);
  int dims = params.getOptionIntValue("-d", 3);
  size_t N = params.getOptionLongValue("-n", -1);
  int tag = params.getOptionIntValue("-t", 1);
  int rounds = params.getOptionIntValue("-r", 3);
  int query_type = params.getOptionIntValue("-q", 0);
  int read_insert_file = params.getOptionIntValue("-i", 1);
  int summary = params.getOptionIntValue("-s", 0);
  int tree_type = params.getOptionIntValue("-T", 0);
  int split_type = params.getOptionIntValue("-l", 0);

  auto test_func = []<class TreeDesc, typename Point>(
                       int const& kDim, parlay::sequence<Point> const& wp,
                       parlay::sequence<Point> const& wi, size_t const& N,
                       int const& K, int const& kRounds,
                       string const& kInsertFile, int const& kTag,
                       int const& kQueryType, int const kSummary) {
    using Tree = TreeDesc::TreeType;
    using Points = typename Tree::Points;

    Tree tree;
    constexpr bool kTestTime = true;

    BuildTree<Point, Tree, kTestTime, 2>(wp, kRounds, tree);

    // NOTE: batch insert
    if (kTag & (1 << 0)) {
      if (kSummary) {
        parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
        for (size_t i = 0; i < ratios.size(); i++) {
          BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, ratios[i]);
        }
      } else {
        BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, 0.1);
      }
    }

    // NOTE: batch delete
    if (kTag & (1 << 1)) {
      if (kSummary) {
        parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
        for (size_t i = 0; i < ratios.size(); i++) {
          BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds, ratios[i]);
        }
      } else {
        BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds,
                                            kBatchInsertRatio);
      }
    }
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
    //     BatchDiff<Point, Tree, kTestTime>(
    //         tree, wp, kRounds, kBatchDiffTotalRatio, kBatchDiffOverlapRatio);
    //   }
    // }
    //
    // // WARN: compress the kdnode to MultiNode, should remove except for
    // // exp if constexpr (IsKdTree<Tree>) {
    // //   tree.Compress2Multi();
    // // }
    //

    if (kTag & (1 << 3)) {
      parlay::sequence<double> const ratios = {0.1, 0.01, 0.001};
      for (auto rat : ratios) {
        BatchUpdateByStep<Point, Tree, true>(tree, wp, wi, kRounds, rat);
      }
    }

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

    if (kQueryType & (1 << 1)) {  // NOTE: range count
      if (!kSummary) {
        int recNum = kRangeQueryNum;
        kdknn = new Typename[recNum];

        // std::cout << std::endl;
        for (int i = 0; i < 3; i++) {
          rangeCountFix<Point>(wp, tree, kdknn, kRounds, i, recNum, kDim);
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

    std::cout << "\n" << std::flush;

    tree.DeleteTree();

    return;
  };

  Wrapper::ApplySpacialFillingCurve(tree_type, dims, split_type, params,
                                    test_func);

  return 0;
}
