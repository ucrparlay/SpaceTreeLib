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
    using Leaf = typename Tree::Leaf;
    using Interior = typename Tree::Interior;

    Tree tree;
    constexpr bool kTestTime = true;

    Typename* kdknn = nullptr;
    auto run_batch_knn = [&](Points const& query_pts, int kth) {
      kdknn = new Typename[query_pts.size()];
      queryKNN<Point>(kDim, query_pts, kRounds, tree, kdknn, kth, true);
      delete[] kdknn;
    };

    // NOTE: batch insert by step
    if (kTag & (1 << 3)) {
      parlay::sequence<double> const ratios = {1, 0.1, 0.01, 0.001, 0.0001};
      // parlay::sequence<double> const ratios = {1, 0.001};
      for (auto rat : ratios) {
        BatchInsertByStep<Point, Tree, true>(tree, wp, kRounds, rat);

        // test knn
        if (static_cast<int>(rat) == 1) continue;

        auto new_root = tree.template MergeUp<Leaf, Interior>(tree.GetRoot());
        tree.SetRoot(new_root);

        int k[3] = {1, 10, 100};

        std::cout << "in-dis-skewed knn time: ";
        size_t batch_size = static_cast<size_t>(wp.size() * kBatchQueryRatio);
        for (int i = 0; i < 3; i++) {
          run_batch_knn(wp.subseq(0, batch_size), k[i]);
        }
        puts("");

        std ::cout << "out-dis-skewed knn time: ";
        for (int i = 0; i < 3; i++) {
          run_batch_knn(wp.subseq(wp.size() - batch_size, wp.size()), k[i]);
        }
        puts("");

        // NOTE: sample points within the whole input datasets
        auto query_pts =
            parlay::pack(wp, parlay::tabulate(wp.size(), [&](size_t i) -> bool {
                           return i % (wp.size() / (batch_size * 2)) == 0;
                         }));

        std::cout << "in-dis-uniform knn time: ";
        for (int i = 0; i < 3; i++) {
          run_batch_knn(parlay::random_shuffle(query_pts.subseq(0, batch_size)),
                        k[i]);
        }
        puts("");

        std::cout << "out-dis-uniform knn time: ";
        for (int i = 0; i < 3; i++) {
          run_batch_knn(parlay::random_shuffle(
                            query_pts.subseq(batch_size, query_pts.size())),
                        k[i]);
        }
        puts("");
      }
    }
  };

  // Wrapper::ApplyOrthogonal(tree_type, dims, split_type, params,
  //                          Wrapper::default_test_func);
  Wrapper::ApplyOrthogonal(tree_type, dims, split_type, params, test_func);

  return 0;
}
