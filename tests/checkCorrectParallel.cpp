#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "cpdd/orth_tree.h"
#include "testFramework.h"
#include <CGAL/Orthtree.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <bits/stdc++.h>
#include <iterator>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <CGAL/Fuzzy_sphere.h>
// using point = PointID<coord, 5>;
using point = PointType<Coord, 2>;
using points = parlay::sequence<point>;

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
int Dim, K, tag, rounds;
bool insert;
int queryType, tree_type;
size_t N;

size_t maxReduceSize = 0;
const size_t batchQuerySize = 1000000;
const double batchInsertCheckRatio = 0.1;

void runCGAL(points& wp, points& wi, Typename* cgknn, int queryNum,
             parlay::sequence<Point_d>& Out) {
    //* cgal
    std::vector<Point_d> _points(N);
    parlay::parallel_for(
        0, N,
        [&](size_t i) {
            _points[i] = Point_d(Dim, std::begin(wp[i].pnt),
                                 (std::begin(wp[i].pnt) + Dim));
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
            _points[j] = Point_d(Dim, std::begin(wi[j].pnt),
                                 (std::begin(wi[j].pnt) + Dim));
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

    if (queryType == 0) {  //* NN
        // size_t S = wp.size();
        size_t S = batchQuerySize;
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, S),
                          [&](const tbb::blocked_range<std::size_t>& r) {
                              for (std::size_t s = r.begin(); s != r.end();
                                   ++s) {
                                  // Neighbor search can be instantiated from
                                  // several threads at the same time
                                  Point_d query(Dim, std::begin(wp[s].pnt),
                                                std::begin(wp[s].pnt) + Dim);
                                  Neighbor_search search(tree, query, K);
                                  Neighbor_search::iterator it = search.end();
                                  it--;
                                  cgknn[s] = it->second;
                              }
                          });
    } else if (queryType == 1) {  //* range count
        size_t n = wp.size();
        parlay::parallel_for(0, queryNum, [&](size_t i) {
            std::vector<Point_d> _ans(n);
            Point_d a(Dim, std::begin(wp[i].pnt), std::end(wp[i].pnt)),
                b(Dim, std::begin(wp[(i + n / 2) % n].pnt),
                  std::end(wp[(i + n / 2) % n].pnt));
            Fuzzy_iso_box fib(a, b, 0.0);
            // auto d = cpdd::ParallelKDtree<point>::p2p_distance( wp[i], wp[( i
            // + n / 2 ) % n],
            //                                                     wp[i].get_dim()
            //                                                     );
            // d = static_cast<coord>( std::sqrt( d ) );
            // if ( i == 0 ) {
            //   LOG << wp[i] << d << ENDL;
            // }
            // Fuzzy_circle fib( a, d );
            size_t cnt = 0;
            counter_iterator<size_t> cnt_iter(cnt);

            auto it = tree.search(cnt_iter, fib);
            // auto it = tree.search( _ans.begin(), fib );
            cgknn[i] = cnt;
            // cgknn[i] = std::distance( _ans.begin(), it );
        });
    } else if (queryType == 2) {  //* range query
        size_t n = wp.size();
        assert(maxReduceSize != 0);
        Out.resize(queryNum * maxReduceSize);
        parlay::parallel_for(0, queryNum, [&](size_t i) {
            Point_d a(Dim, std::begin(wp[i].pnt), std::end(wp[i].pnt)),
                b(Dim, std::begin(wp[(i + n / 2) % n].pnt),
                  std::end(wp[(i + n / 2) % n].pnt));
            Fuzzy_iso_box fib(a, b, 0.0);
            auto it = tree.search(Out.begin() + i * maxReduceSize, fib);
            cgknn[i] = std::distance(Out.begin() + i * maxReduceSize, it);
        });
    }

    if (tag == 1) {
        wp.pop_tail(wi.size() * batchInsertCheckRatio);
        assert(wp.size() == N);
    }
    tree.clear();
}

template<typename Tree>
void runKDParallel(points& wp, const points& wi, Typename* kdknn, points& p,
                   int queryNum) {
    puts("build kd tree");
    using pkdtree = Tree;
    using box = typename Tree::Box;
    pkdtree pkd;
    size_t n = wp.size();

    buildTree<point>(Dim, wp, rounds, pkd);
    cpdd::Node* KDParallelRoot = pkd.GetRoot();
    pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                          typename Tree::SplitRuleType>();

    if (tag >= 1) {
        BatchInsert<point, Tree>(pkd, wp, wi, Dim, 2, batchInsertCheckRatio);
        if (tag == 1) wp.append(wi.cut(0, wp.size() * batchInsertCheckRatio));
        pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                              typename Tree::SplitRuleType>();
        LOG << "finish insert" << ENDL;
    }

    if (tag >= 2) {
        batchDelete<point, Tree, kBatchDiff>(pkd, wp, wi, Dim, 2, true,
                                             batchInsertCheckRatio);
        pkd.template Validate<typename Tree::Leaf, typename Tree::Interior,
                              typename Tree::SplitRuleType>();
        LOG << "finish delete" << ENDL;
    }

    //* query phase

    assert(N >= K);
    assert(tag == 1 || wp.size() == N);
    LOG << "begin kd query" << ENDL;
    if (queryType == 0) {
        points new_wp(batchQuerySize);
        parlay::copy(wp.cut(0, batchQuerySize), new_wp.cut(0, batchQuerySize));
        queryKNN<point>(Dim, wp, rounds, pkd, kdknn, K, true);
    } else if (queryType == 1) {
        rangeCount<point>(wp, pkd, kdknn, rounds, queryNum);
        // rangeCountRadius<point>( wp, pkd, kdknn, rounds, queryNum );
    } else if (queryType == 2) {
        rangeCount<point>(wp, pkd, kdknn, rounds, queryNum);
        maxReduceSize =
            parlay::reduce(parlay::delayed_tabulate(
                               queryNum, [&](size_t i) { return kdknn[i]; }),
                           parlay::maximum<Typename>());
        LOG << "max size " << maxReduceSize << ENDL;
        p.resize(queryNum * maxReduceSize);
        rangeQuery<point>(wp, pkd, kdknn, rounds, queryNum, p);
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
    char* iFile = P.getOptionValue("-p");
    K = P.getOptionIntValue("-k", 100);
    Dim = P.getOptionIntValue("-d", 3);
    N = P.getOptionLongValue("-n", -1);
    tag = P.getOptionIntValue("-t", 0);
    rounds = P.getOptionIntValue("-r", 3);
    queryType = P.getOptionIntValue("-q", 0);
    char* _insertFile = P.getOptionValue("-i");
    tree_type = P.getOptionIntValue("-T", 1);

    assert(Dim == point().get_dim());

    if (tag == 0) {
        puts("build and query");
    } else if (tag == 1) {
        puts("build, insert and query");
    }

    int LEAVE_WRAP = 32;
    points wp;
    std::string name, insertFile;

    //* initialize points
    if (iFile != NULL) {
        name = std::string(iFile);
        name = name.substr(name.rfind("/") + 1);
        std::cout << name << " ";
        auto [n, d] = read_points<point>(iFile, wp, K);
        N = n;
        assert(Dim == d);
    } else {  //* construct data byself
        K = 100;
        generate_random_points<point>(wp, 10000, N, Dim);
        assert(wp.size() == N);
        std::string name =
            std::to_string(N) + "_" + std::to_string(Dim) + ".in";
        std::cout << name << " ";
    }

    LOG << std::setprecision(13) << wp[0] << wp[1] << ENDL;

    using ref_t = std::reference_wrapper<point>;
    using pairs = std::pair<ref_t, int>;
    ref_t a(wp[0]);
    pairs p(a, 0);
    LOG << p.first << " " << p.second << ENDL;

    Typename* cgknn;
    Typename* kdknn;
    points wi;
    int queryNum = N / 10000;

    //* initialize insert points file
    if (tag >= 1 && iFile != NULL) {
        if (_insertFile == NULL) {
            int id = std::stoi(name.substr(0, name.find_first_of('.')));
            id = (id + 1) % 10;  //! MOD graph number used to test
            if (!id) id++;
            int pos = std::string(iFile).rfind("/") + 1;
            insertFile =
                std::string(iFile).substr(0, pos) + std::to_string(id) + ".in";
        } else {
            insertFile = std::string(_insertFile);
        }
        std::cout << insertFile << ENDL;
    }

    //* generate points
    if (tag >= 1) {
        if (iFile == NULL) {
            generate_random_points<point>(wi, 1000000, N / 2, Dim);
            LOG << "insert " << N / 5 << " points" << ENDL;
        } else {
            auto [nn, nd] = read_points<point>(insertFile.c_str(), wi, K);
            if (nd != Dim || nn != N) {
                puts("read inserted points dimension wrong");
                abort();
            } else {
                puts("read inserted points from file");
            }
        }
    }

    //* set result array size
    if (queryType == 0) {  //*NN
        LOG << "---do NN query---" << ENDL;
        if (tag == 0) {
            cgknn = new Typename[N];
            kdknn = new Typename[N];
        } else if (tag == 1) {
            puts("insert points from file");
            cgknn = new Typename[N + wi.size()];
            kdknn = new Typename[N + wi.size()];
        } else if (tag == 2) {
            puts("insert then delete points from file");
            cgknn = new Typename[N];
            kdknn = new Typename[N];
        }
    } else if (queryType == 1) {  //* range Count
        LOG << "---do range Count---" << ENDL;

        cgknn = new Typename[queryNum];
        kdknn = new Typename[queryNum];
    } else if (queryType == 2) {
        LOG << "---do range Query---" << ENDL;

        cgknn = new Typename[queryNum];
        kdknn = new Typename[queryNum];
    }

    points kdOut;
    parlay::sequence<Point_d> cgOut;
    if (tree_type == 0) {
        LOG << "test kd tree" << ENDL;
        runKDParallel<cpdd::KdTree<point, cpdd::MaxStretchDim<point>>>(
            // runKDParallel<cpdd::KdTree<point, cpdd::RotateDim<point>>>(
            wp, wi, kdknn, kdOut, queryNum);
    }
    // else if (tree_type == 1 && Dim == 2) {
    //     LOG << "test quad tree" << ENDL;
    //     runKDParallel<cpdd::OrthTree<point, cpdd::RotateDim<point>, 2>>(
    //         wp, wi, kdknn, kdOut, queryNum);
    // } else if (tree_type == 1 && Dim == 3) {
    //     LOG << "test oct tree" << ENDL;
    //     runKDParallel<cpdd::OrthTree<point, cpdd::RotateDim<point>, 3>>(
    //         wp, wi, kdknn, kdOut, queryNum);
    // }
    runCGAL(wp, wi, cgknn, queryNum, cgOut);

    //* verify
    if (queryType == 0) {
        LOG << "check NN" << ENDL;
        size_t S = batchQuerySize;
        for (size_t i = 0; i < S; i++) {
            if (std::abs(cgknn[i] - kdknn[i]) > 1e-4) {
                puts("");
                puts("wrong");
                std::cout << i << " " << cgknn[i] << " " << kdknn[i]
                          << std::endl;
                return 0;
            }
        }
    } else if (queryType == 1) {
        LOG << "check range count" << ENDL;
        for (size_t i = 0; i < queryNum; i++) {
            if (std::abs(cgknn[i] - kdknn[i]) > 1e-4) {
                puts("");
                puts("wrong");
                std::cout << i << " " << cgknn[i] << " " << kdknn[i]
                          << std::endl;
                //        return 0;
            }
        }
    } else if (queryType == 2) {
        LOG << "check range query" << ENDL;
        assert(kdOut.size() == cgOut.size());
        auto kdans = parlay::tabulate(kdOut.size(), [&](size_t i) {
            return Point_d(Dim, std::begin(kdOut[i].pnt),
                           (std::begin(kdOut[i].pnt) + Dim));
        });

        // parlay::parallel_for( 0, queryNum, [&]( size_t i ) {
        for (int i = 0; i < queryNum; i++) {
            if (std::abs(cgknn[i] - kdknn[i]) > 1e-4) {
                puts("");
                puts("count num wrong");
                std::cout << i << " " << cgknn[i] << " " << kdknn[i]
                          << std::endl;
                return 0;
            }

            size_t s = i * maxReduceSize;
            parlay::sort_inplace(kdans.cut(s, s + size_t(kdknn[i])));
            parlay::sort_inplace(cgOut.cut(s, s + size_t(cgknn[i])));
            // std::sort(kdans.begin() + s, kdans.begin() + s +
            // size_t(kdknn[i])); std::sort(cgOut.begin() + s, cgOut.begin() + s
            // + size_t(cgknn[i]));

            for (int j = 0; j < kdknn[i]; j++) {
                if (kdans[j + s] != cgOut[j + s]) {
                    puts("");
                    puts("point wrong");
                    std::cout << j << " " << cgOut[j + s] << " " << kdans[j + s]
                              << std::endl;
                    return 0;
                }
            }
        }
        // } );
    }

    puts("\nok");
    return 0;
}
