#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "test_framework.h"

template <typename Point, typename Tree, bool printHeight = 0,
          bool printVisNode = 1>
void QueryFind(Tree& tree, uint_fast8_t const& Dim,
               parlay::sequence<Point> const& WP, int const& rounds,
               Typename* kdknn, int const K) {
  using Points = typename Tree::Points;
  using Coord = typename Point::Coord;
  using DisType = typename Point::DisType;
  using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  using SplitRule = typename Tree::SplitRule;
  size_t n = WP.size();

  // Points wp = WP;
  Points wp = parlay::random_shuffle(WP);
  auto wp_id = parlay::tabulate(n, [&](size_t i) {
    return std::make_pair(SplitRule::Encode(wp[i]), wp[i].GetAug().id);
  });

  double aveQuery = time_loop(
      rounds, 0.1, [&]() {},
      [&]() {
        parlay::parallel_for(0, n, [&](size_t i) { tree.find(wp_id[i]); });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;

  return;
}
template <typename Point, typename Tree>
void OneDimensionalFind(Tree& tree, uint_fast8_t const& Dim,
                        parlay::sequence<Point> const& WP, int const& rounds) {
  using Points = typename Tree::Points;
  using Coord = typename Point::Coord;
  using DisType = typename Point::DisType;
  using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  using SplitRule = typename Tree::SplitRule;
  size_t n = WP.size();

  Points wp = WP;
  // Points wp = parlay::random_shuffle(WP);
  auto wp_id = parlay::tabulate(n, [&](size_t i) { return wp[i].aug.id; });

  double aveQuery = time_loop(
      rounds, 0.1, [&]() {},
      [&]() {
        parlay::parallel_for(0, n, [&](size_t i) { tree.find(wp_id[i]); });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;

  return;
}

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

  auto test_func = []<class TreeDesc, typename Point>(
                       int const& kDim, parlay::sequence<Point> const& wp,
                       parlay::sequence<Point> const& wi, size_t const& N,
                       int const& K, int const& kRounds,
                       string const& kInsertFile, int const& kTag,
                       int const& kQueryType, int const kSummary) {
    using Tree = TreeDesc::TreeType;
    using Points = typename Tree::Points;
    using Box = typename Tree::Box;
    Tree tree;
    constexpr bool kTestTime = true;
    auto new_wp =
        parlay::random_shuffle(parlay::tabulate(wp.size(), [&](size_t i) {
          Point pt = wp[i];
          pt.pnt[0] = i;
          return pt;
        }));

    // std::cout << "begin build tree: ";
    BuildTree<Point, Tree, kTestTime, 2>(new_wp, kRounds, tree);
    // std::cout << "begin query find: ";
    // QueryFind<Point, Tree, 0, 1>(tree, kDim, wp, kRounds,
    //                              static_cast<Typename*>(nullptr), K);
    OneDimensionalFind<Point, Tree>(tree, kDim, new_wp, kRounds);
  };

  Wrapper::ApplyBaselines(tree_type, dims, split_type, params, test_func);

  return 0;
}
