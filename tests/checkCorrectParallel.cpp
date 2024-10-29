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

#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "cpdd/orth_tree.h"
#include "testFramework.h"
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
int query_num = 10000;
size_t const batchQuerySize = 1000000;
double const batchInsertCheckRatio = 0.1;

template <typename Point>
void runCGAL(auto& wp, auto& wi, Typename* cgknn,
             [[maybe_unused]] int query_num,
             [[maybe_unused]] parlay::sequence<Point_d>& Out, int const tag,
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

  // LOG << tree.bounding_box() << ENDL;
  size_t sz = wp.size() * batchInsertCheckRatio;

  if (tag >= 1) {
    _points.resize(wi.size());
    parlay::parallel_for(0, wi.size(), [&](size_t j) {
      _points[j] =
          Point_d(kDim, std::begin(wi[j].pnt), (std::begin(wi[j].pnt) + kDim));
    });
    tree.insert(_points.begin(), _points.begin() + sz);
    tree.build<CGAL::Parallel_tag>();
    assert(tree.size() == wp.size() + sz);
    wp.append(wi.cut(0, sz));
    puts("finish insert to cgal");
  }
  LOG << tree.root()->num_items() << ENDL;

  if (tag >= 2) {
    assert(_points.size() == wi.size());
    // for (auto p : _points) {
    for (size_t i = 0; i < sz; i++) {
      tree.remove(_points[i]);
    }

    LOG << tree.root()->num_items() << ENDL;
    wp.pop_tail(sz);
    puts("finish delete from cgal");
  }

  //* cgal query
  LOG << "begin tbb query" << ENDL << std::flush;
  assert(tree.is_built());

  if (query_type == 0) {  //* NN
    // size_t S = wp.size();
    size_t S = batchQuerySize;
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
  } else if (query_type == 1) {  //* range count
    // size_t n = wp.size();
    // parlay::parallel_for(0, query_num, [&](size_t i) {
    //   std::vector<Point_d> _ans(n);
    //   Point_d a(kDim, std::begin(wp[i].pnt), std::end(wp[i].pnt)),
    //       b(kDim, std::begin(wp[(i + n / 2) % n].pnt),
    //         std::end(wp[(i + n / 2) % n].pnt));
    //   Fuzzy_iso_box fib(a, b, 0.0);
    //   // auto d = cpdd::ParallelKDtree<point>::p2p_distance( wp[i], wp[( i
    //   // + n / 2 ) % n],
    //   //                                                     wp[i].GetDim()
    //   //                                                     );
    //   // d = static_cast<coord>( std::sqrt( d ) );
    //   // if ( i == 0 ) {
    //   //   LOG << wp[i] << d << ENDL;
    //   // }
    //   // Fuzzy_circle fib( a, d );
    //   size_t cnt = 0;
    //   counter_iterator<size_t> cnt_iter(cnt);
    //
    //   [[maybe_unused]] auto it = tree.search(cnt_iter, fib);
    //   // auto it = tree.search( _ans.begin(), fib );
    //   cgknn[i] = cnt;
    // cgknn[i] = std::distance( _ans.begin(), it );
    // });
  } else if (query_type == 2) {  //* range query
    // size_t n = wp.size();
    // assert(maxReduceSize != 0);
    // Out.resize(query_num * maxReduceSize);
    // parlay::parallel_for(0, query_num, [&](size_t i) {
    //   Point_d a(kDim, std::begin(wp[i].pnt), std::end(wp[i].pnt)),
    //       b(kDim, std::begin(wp[(i + n / 2) % n].pnt),
    //         std::end(wp[(i + n / 2) % n].pnt));
    //   Fuzzy_iso_box fib(a, b, 0.0);
    //   auto it = tree.search(Out.begin() + i * maxReduceSize, fib);
    //   cgknn[i] = std::distance(Out.begin() + i * maxReduceSize, it);
    // });
  }

  if (tag == 1) {
    wp.pop_tail(wi.size() * batchInsertCheckRatio);
    assert(wp.size() == N);
  }
  tree.clear();
}

template <typename Point, typename TreeDesc>
void runKDParallel(auto& wp, auto const& wi, Typename* kdknn, auto& p,
                   int query_num, int const query_type, int const K,
                   int const tag, int const rounds) {
  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;
  size_t const N = wp.size();
  size_t const kDim = std::tuple_size_v<typename Point::Coords>;

  Tree pkd;
  // [[maybe_unused]] size_t n = wp.size();

  puts("build tree");
  buildTree<Point>(kDim, wp, rounds, pkd);
  // [[maybe_unused]] cpdd::Node* KDParallelRoot = pkd.GetRoot();
  pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                        typename Tree::SplitRuleType>();

  if (tag >= 1) {
    BatchInsert<Point, Tree>(pkd, wp, wi, kDim, 2, batchInsertCheckRatio);
    if (tag == 1) wp.append(wi.cut(0, wp.size() * batchInsertCheckRatio));
    pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                          typename Tree::SplitRuleType>();
    LOG << "finish insert" << ENDL;
  }

  if (tag >= 2) {
    batchDelete<Point, Tree, kBatchDelete>(pkd, wp, wi, kDim, 2, true,
                                           batchInsertCheckRatio);
    pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                          typename Tree::SplitRuleType>();
    LOG << "finish delete" << ENDL;
  }

  // NOTE: query phase
  assert(tag == 1 || wp.size() == N);
  LOG << "begin kd query" << ENDL;
  if (query_type == 0) {
    Points new_wp(batchQuerySize);
    parlay::copy(wp.cut(0, batchQuerySize), new_wp.cut(0, batchQuerySize));
    queryKNN<Point>(kDim, wp, rounds, pkd, kdknn, K, true);
  } else if (query_type == 1) {
    for (auto range_query_type : {0, 1, 2}) {
      rangeCount<Point>(wp, pkd, kdknn, rounds, query_num, range_query_type,
                        kDim);
    }
  } else if (query_type == 2) {
    for (auto range_query_type : {0, 1, 2}) {
      rangeQuery<Point>(wp, pkd, kdknn, rounds, query_num, range_query_type,
                        kDim, p);
    }
  }

  if (tag == 1) wp.pop_tail(wi.size() * batchInsertCheckRatio);
  assert(wp.size() == N);
  pkd.DeleteTree();
  return;
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>]");
  char* input_file_path = P.getOptionValue("-p");
  int K = P.getOptionIntValue("-k", 100);
  int kDim = P.getOptionIntValue("-d", 3);
  size_t N = P.getOptionLongValue("-n", -1);
  int tag = P.getOptionIntValue("-t", 1);
  int kRounds = P.getOptionIntValue("-r", 3);
  int query_type = P.getOptionIntValue("-q", 0);
  int read_insert_file = P.getOptionIntValue("-i", 1);
  int tree_type = P.getOptionIntValue("-T", 0);

  if (tag == 0) {
    puts("build and query");
  } else if (tag == 1) {
    puts("build, insert and query");
  }

  auto run_test = [&]<class Wrapper>(Wrapper) {
    auto run = [&](auto dim_wrapper) -> void {
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
        id = (id + 1) % 10;  // WARN: MOD graph number used to test
        if (!id) id++;
        int pos = std::string(input_file_path).rfind("/") + 1;
        insert_file_path = std::string(input_file_path).substr(0, pos) +
                           std::to_string(id) + ".in";
        auto [n, d] =
            read_points<PointTypeAlias>(insert_file_path.c_str(), wi, K);
        assert(d == kDim);
      }

      // NOTE: begin the test
      // NOTE: alloc the memory
      Coord* cgknn;
      Coord* kdknn;
      if (query_type == 0) {  //*NN
        LOG << "---do NN query---" << ENDL;
        if (tag == 0) {
          cgknn = new Coord[wp.size()];
          kdknn = new Coord[wp.size()];
        } else if (tag == 1) {
          puts("insert points from file");
          cgknn = new Coord[wp.size() + wi.size()];
          kdknn = new Coord[wp.size() + wi.size()];
        } else if (tag == 2) {
          puts("insert then delete points from file");
          cgknn = new Coord[wp.size()];
          kdknn = new Coord[wp.size()];
        } else {
          puts("wrong tag");
          abort();
        }
      } else if (query_type == 1) {  //* range Count
        LOG << "---do range Count---" << ENDL;
        kdknn = new Coord[query_num];
      } else if (query_type == 2) {
        LOG << "---do range Query---" << ENDL;
        kdknn = new Coord[query_num];
      } else {
        puts("wrong query type");
        abort();
      }

      // NOTE: run the test
      Points kd_range_out;
      parlay::sequence<Point_d> cg_range_out;

      runKDParallel<PointTypeAlias, Desc>(
          wp, wi, kdknn, kd_range_out, query_num, query_type, K, tag, kRounds);

      if (query_type == 0) {
        runCGAL<PointTypeAlias>(wp, wi, cgknn, query_num, cg_range_out, tag,
                                query_type, K);
      }

      // NOTE: verify
      if (query_type == 0) {
        LOG << "check NN" << ENDL;
        size_t S = batchQuerySize;
        for (size_t i = 0; i < S; i++) {
          if (std::abs(cgknn[i] - kdknn[i]) > 1e-4) {
            puts("");
            puts("wrong");
            std::cout << i << " " << cgknn[i] << " " << kdknn[i] << std::endl;
            abort();
          }
        }
      }
    };

    if (tag == -1) {
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
    LOG << "run KDtree" << ENDL;
    run_test(wrapper::KDtree{});
  } else if (tree_type == 1 && kDim == 2) {
    LOG << "run QuadTree" << ENDL;
    run_test(wrapper::QuadTree{});
  } else if (tree_type == 1 && kDim == 3) {
    LOG << "run OctTree" << ENDL;
    run_test(wrapper::OctTree{});
  }

  puts("\nok");
  return 0;
}
