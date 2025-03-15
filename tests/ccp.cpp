#include <CGAL/Cartesian_d.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_d.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <iterator>
#include <tuple>

#include "pspt/dependence/splitter.h"
#include "pspt/kd_tree.h"
#include "pspt/orth_tree.h"
#include "test_framework.h"
// using point = PointID<coord, 5>;
// using point = PointType<Coord, 2>;
// using points = parlay::sequence<point>;

typedef CGAL::Cartesian_d<Typename> Kernel;
typedef Kernel::Point_d Point_d;
typedef CGAL::Search_traits_d<Kernel> TreeTraits;
typedef CGAL::Random_points_in_cube_d<Point_d> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>
    N_Random_points_iterator;
typedef CGAL::Median_of_rectangle<TreeTraits> Median_of_rectangle;
typedef CGAL::Euclidean_distance<TreeTraits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits, Distance,
                                           Median_of_rectangle>
    Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef CGAL::Fuzzy_iso_box<TreeTraits> Fuzzy_iso_box;
typedef CGAL::Fuzzy_sphere<TreeTraits> Fuzzy_circle;

size_t maxReduceSize = 0;
int const kCCPQueryNum = 10000;
size_t const kCCPBatchQuerySize = 1000000;
// size_t const kCCPBatchQuerySize = 110;
double const batchInsertCheckRatio = 0.1;
double const kCCPBatchDiffTotalRatio = 1.0;
double const kCCPBatchDiffOverlapRatio = 0.5;

template <typename Point>
void runCGAL(auto const& wp, auto const& wi, Typename* cgknn,
             [[maybe_unused]] int query_num, int const tag,
             int const query_type, int const K) {
  //* cgal
  size_t const N = wp.size();
  size_t const kDim = std::tuple_size_v<typename Point::Coords>;
  std::vector<Point_d> _points(N);
  parlay::parallel_for(
      0, N,
      [&](size_t i) {
        _points[i] = Point_d(kDim, std::begin(wp[i].pnt),
                             (std::begin(wp[i].pnt) + kDim));
      },
      1000);
  Median_of_rectangle median;
  Tree tree(_points.begin(), _points.end(), median);
  tree.build<CGAL::Parallel_tag>();
  // std::cout << tree.bounding_box() << std::endl;
  size_t sz = wp.size() * batchInsertCheckRatio;

  if (tag == 1) {  // NOTE: only insert
    _points.resize(wi.size());
    parlay::parallel_for(0, wi.size(), [&](size_t j) {
      _points[j] =
          Point_d(kDim, std::begin(wi[j].pnt), (std::begin(wi[j].pnt) + kDim));
    });
    tree.insert(_points.begin(), _points.begin() + sz);
    tree.build<CGAL::Parallel_tag>();
    assert(tree.size() == wp.size() + sz);
    puts("finish insert to cgal");
  }
  std::cout << tree.root()->num_items() << std::endl;

  if (tag & (1 << 2)) {
    size_t total_batch_size = static_cast<size_t>(N * kCCPBatchDiffTotalRatio);
    size_t overlap_size =
        static_cast<size_t>(total_batch_size * kCCPBatchDiffOverlapRatio);
    for (size_t i = 0; i < overlap_size; i++) {
      tree.remove(_points[i]);
    }
    puts("finish diff from cgal");
  }

  //* cgal query
  std::cout << "begin tbb query" << std::endl << std::flush;
  assert(tree.is_built());

  if (query_type & (1 << 0)) {  //* NN
    // size_t S = wp.size();
    size_t S = kCCPBatchQuerySize;
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, S),
                      [&](tbb::blocked_range<std::size_t> const& r) {
                        for (std::size_t s = r.begin(); s != r.end(); ++s) {
                          // Neighbor search can be instantiated from
                          // several threads at the same time
                          Point_d query(kDim, std::begin(wp[s].pnt),
                                        std::begin(wp[s].pnt) + kDim);
                          Neighbor_search search(tree, query, K);
                          Neighbor_search::iterator it = search.end();
                          it--;
                          cgknn[s] = it->second;
                        }
                      });
  }

  tree.clear();
}

template <typename Point, typename TreeDesc>
void runKDParallel(auto const& wp, auto const& wi, Typename* kdknn,
                   Typename* cgknn, int query_num, int const query_type,
                   int const K, int const tag, int const rounds) {
  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;
  size_t const kDim = std::tuple_size_v<typename Point::Coords>;

  Tree tree;
  constexpr bool kTestTime = false;

  BuildTree<Point, Tree, kTestTime>(wp, rounds, tree);
  tree.template Validate<typename Tree::Leaf, typename Tree::Interior,
                         typename Tree::SplitRuleType>();

  puts("---------------finish build tree---------------");

  // if (tag & (1 << 0)) {
  //   BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, 2,
  //   batchInsertCheckRatio); tree.template Validate<typename Tree::Leaf,
  //   typename Tree::Interior,
  //                          typename Tree::SplitRuleType>();
  //   std::cout << "---------------finish insert----------------\n" <<
  //   std::flush;
  // }
  //
  // if (tag &
  //     (1 << 1)) {  // NOTE: without tag & (1<<0), this will remove no points
  //   BatchDelete<Point, Tree, kTestTime>(tree, wp, wi, 2,
  //   batchInsertCheckRatio); tree.template Validate<typename Tree::Leaf,
  //   typename Tree::Interior,
  //                          typename Tree::SplitRuleType>();
  //   std::cout << "---------------finish delete----------------\n" <<
  //   std::flush;
  // }
  //
  // if (tag & (1 << 2)) {
  //   BatchDiff<Point, Tree, kTestTime>(tree, wp, 2, kCCPBatchDiffTotalRatio,
  //                                     kCCPBatchDiffOverlapRatio);
  //   // BatchDiff<Point, Tree, kTestTime>(tree, wp, 2, 0.5, 1.0);
  //   tree.template Validate<typename Tree::Leaf, typename Tree::Interior,
  //                          typename Tree::SplitRuleType>();
  //   assert(tree.GetRoot()->size ==
  //          wp.size() - static_cast<size_t>(wp.size() *
  //          kCCPBatchDiffTotalRatio *
  //                                          kCCPBatchDiffOverlapRatio));
  //   std::cout << "---------------finish diff------------------\n" <<
  //   std::flush;
  // }

  // NOTE: query phase
  if (query_type & (1 << 0)) {  // NOTE: NN query
    Points new_wp(kCCPBatchQuerySize);
    parlay::copy(wp.cut(0, kCCPBatchQuerySize),
                 new_wp.cut(0, kCCPBatchQuerySize));
    queryKNN<Point>(kDim, new_wp, rounds, tree, kdknn, K, true);
    std::cout << "run cgal\n" << std::flush;
    runCGAL<Point>(wp, wi, cgknn, query_num, tag, query_type, K);
    std::cout << "check NN\n" << std::flush;
    size_t S = kCCPBatchQuerySize;
    for (size_t i = 0; i < S; i++) {
      if (std::abs(cgknn[i] - kdknn[i]) > 1e-4) {
        puts("");
        puts("wrong");
        std::cout << i << " " << cgknn[i] << " " << kdknn[i] << "\n"
                  << std::flush;
        abort();
      }
    }
    std::cout << "--------------finish NN query------------------\n"
              << std::flush;
  }

  // if (query_type & (1 << 1)) {  // NOTE: range count
  //   for (auto range_query_type : {0, 1, 2}) {
  //     if (tag & (1 << 2)) {
  //       Points new_wp(tree.GetRoot()->size);
  //       tree.Flatten(new_wp);
  //       RangeCount<Point>(new_wp, tree, rounds, query_num, range_query_type,
  //                         kDim);
  //     } else {
  //       RangeCount<Point>(wp, tree, rounds, query_num, range_query_type,
  //       kDim);
  //     }
  //   }
  //   std::cout << "--------------finish range count------------------\n"
  //             << std::flush;
  // }
  //
  // if (query_type & (1 << 2)) {  // NOTE: range query
  //   for (auto range_query_type : {0, 1, 2}) {
  //     if (tag & (1 << 2)) {
  //       Points new_wp(tree.GetRoot()->size);
  //       tree.Flatten(new_wp);
  //       RangeQuery<Point>(new_wp, tree, rounds, kCCPQueryNum,
  //       range_query_type,
  //                         kDim);
  //     } else {
  //       RangeQuery<Point>(wp, tree, rounds, kCCPQueryNum, range_query_type,
  //                         kDim);
  //     }
  //   }
  //   std::cout << "--------------finish range query------------------\n"
  //             << std::flush;
  // }

  tree.DeleteTree();
  return;
}

int main(int argc, char* argv[]) {
  commandLine params(argc, argv);

  char* input_file_path = params.getOptionValue("-p");
  int K = params.getOptionIntValue("-k", 100);
  int dim = params.getOptionIntValue("-d", 3);
  size_t N = params.getOptionLongValue("-n", -1);
  int tag = params.getOptionIntValue("-t", 1);
  int kRounds = params.getOptionIntValue("-r", 3);
  int query_type = params.getOptionIntValue("-q", 0);
  int read_insert_file = params.getOptionIntValue("-i", 1);
  int tree_type = params.getOptionIntValue("-T", 0);
  int split_type = params.getOptionIntValue("-l", 0);

  auto run = [&]<typename TreeWrapper, typename Point>(
                 int const& kDim, parlay::sequence<Point> const& wp,
                 parlay::sequence<Point> const& wi, size_t const& N,
                 int const& K, int const& kRounds, string const& kInsertFile,
                 int const& kTag, int const& kQueryType,
                 int const kSummary) -> void {
    // using Point = typename TreeWrapper::Point;
    // using Points = parlay::sequence<Point>;
    // constexpr auto kDim = Point::GetDim();
    //
    // PrintTreeParam<TreeWrapper>();
    //
    // std::string name, insert_file_path = "";
    // Points wp, wi;
    //
    // if (input_file_path != NULL) {  // NOTE: read main Points
    //   name = std::string(input_file_path);
    //   name = name.substr(name.rfind('/') + 1);
    //   std::cout << name << "\n";
    //   [[maybe_unused]] auto [n, d] = read_points<Point>(input_file_path, wp,
    //   K); N = n;
    //   // wp = wp.subseq(0, kCCPBatchQuerySize);
    //   assert(d == kDim);
    // }
    //
    // if (read_insert_file == 1) {  // NOTE: read Points to be inserted
    //   int id = std::stoi(name.substr(0, name.find_first_of('.')));
    //   id = (id + 1) % 10;  // WARN: MOD graph number used to test
    //   if (!id) id++;
    //   auto pos = std::string(input_file_path).rfind('/') + 1;
    //   insert_file_path = std::string(input_file_path).substr(0, pos) +
    //                      std::to_string(id) + ".in";
    //   [[maybe_unused]] auto [n, d] =
    //       read_points<Point>(insert_file_path.c_str(), wi, K);
    //   assert(N == n && d == kDim);
    // }

    // NOTE: begin the test
    // NOTE: alloc the memory
    Coord* cgknn;
    Coord* kdknn;
    if (query_type & (1 << 0)) {  //*NN
      if (tag == 0) {
        cgknn = new Coord[wp.size()];
        kdknn = new Coord[wp.size()];
      } else if (tag & (1 << 0)) {
        cgknn = new Coord[wp.size()];
        kdknn = new Coord[wp.size()];
      } else if (tag & (1 << 1)) {
        cgknn = new Coord[wp.size() + wi.size()];
        kdknn = new Coord[wp.size() + wi.size()];
      } else if (tag & (1 << 2)) {
        cgknn = new Coord[wp.size()];
        kdknn = new Coord[wp.size()];
      } else {
        puts("wrong tag");
        abort();
      }
    } else if (query_type & (1 << 1)) {  //* range Count
      kdknn = new Coord[1];
      cgknn = new Coord[1];
    } else if (query_type & (1 << 2)) {
      kdknn = new Coord[1];
      cgknn = new Coord[1];
    } else {
      puts("Nothing to be queried");
      // abort();
    }

    // NOTE: run the test
    runKDParallel<Point, TreeWrapper>(wp, wi, kdknn, cgknn, kCCPQueryNum,
                                      query_type, K, tag, kRounds);
    delete[] cgknn;
    delete[] kdknn;
  };

  if (tree_type == 0 || tree_type == 1) {
    Wrapper::ApplyOrthogonal(tree_type, dim, split_type, params, run);
  } else if (tree_type == 2) {
    // Wrapper::ApplySpacialFillingCurve(tree_type, dim, split_type, params,
    // run);
  }

  puts("\nok");
  return 0;
}
