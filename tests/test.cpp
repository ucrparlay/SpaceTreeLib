#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "parlay/primitives.h"
#include "testFramework.h"

template<class TreeDesc, typename Point>
void TestSpacialTree(const int& kDim, const parlay::sequence<Point>& wp,
                     const parlay::sequence<Point>& wi, const size_t& N,
                     const int& K, const int& kRounds,
                     const string& kInsertFile, const int& kTag,
                     const int& kQueryType, const int kSummary) {
    using Tree = TreeDesc::TreeType;
    using Points = typename Tree::Points;

    Tree tree;

    buildTree<Point, Tree>(kDim, wp, kRounds, tree);

    Typename* kdknn = nullptr;

    // //* batch insert
    // if (tag >= 1) {
    //     if (summary) {
    //         const parlay::sequence<double> ratios = {0.0001, 0.001, 0.01,
    //         0.1}; for (int i = 0; i < ratios.size(); i++) {
    //             batchInsert<point>(pkd, wp, wi, Dim, rounds, ratios[i]);
    //         }
    //     } else {
    //         batchInsert<point>(pkd, wp, wi, Dim, rounds, batchInsertRatio);
    //     }
    // }
    //
    // //* batch delete
    // if (tag >= 2) {
    //     if (summary) {
    //         const parlay::sequence<double> ratios = {0.0001, 0.001, 0.01,
    //         0.1}; for (int i = 0; i < ratios.size(); i++) {
    //             batchDelete<point>(pkd, wp, wi, Dim, rounds, 0, ratios[i]);
    //         }
    //     } else {
    //         batchDelete<point>(pkd, wp, wi, Dim, rounds, 0,
    //         batchInsertRatio);
    //     }
    // }
    //
    // if (kQueryType & (1 << 0)) {  // NOTE: KNN
    //     auto run_batch_knn = [&](const Points& pts, int kth, size_t
    //     batchSize) {
    //         Points newPts(batchSize);
    //         parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
    //         kdknn = new Typename[batchSize];
    //         queryKNN<Point>(kDim, newPts, kRounds, tree, kdknn, kth, true);
    //         delete[] kdknn;
    //     };
    //
    //     size_t batchSize = static_cast<size_t>(wp.size() * batchQueryRatio);
    //
    //     if (kSummary == 0) {
    //         int k[3] = {1, 10, 100};
    //         for (int i = 0; i < 3; i++) {
    //             run_batch_knn(wp, k[i], batchSize);
    //         }
    //     } else {  // test summary
    //         run_batch_knn(wp, K, batchSize);
    //     }
    // }
    //
    // if (queryType & (1 << 1)) {  // NOTE: batch NN query
    //
    //     // auto run_batch_knn = [&](const points& pts, size_t
    //     // batchSize) {
    //     //   points newPts(batchSize);
    //     //   parlay::copy(pts.cut(0, batchSize), newPts.cut(0,
    //     //   batchSize)); kdknn = new Typename[batchSize];
    //     //   queryKNN<point, true, true>(Dim, newPts, rounds, pkd,
    //     //   kdknn, K, true); delete[] kdknn;
    //     // };
    //     //
    //     // run_batch_knn(wp, static_cast<size_t>(wp.size() *
    //     // batchQueryRatio));
    //     // const std::vector<double> batchRatios = {0.001, 0.01, 0.1, 0.2,
    //     0.5};
    //     // for (auto ratio : batchRatios) {
    //     //   run_batch_knn(wp, static_cast<size_t>(wp.size() * ratio));
    //     // }
    //     // for (auto ratio : batchRatios) {
    //     //   run_batch_knn(wi, static_cast<size_t>(wi.size() * ratio));
    //     // }
    // }
    //
    // if (queryType & (1 << 2)) {  // NOTE: range count
    //     int recNum = rangeQueryNum;
    //     kdknn = new Typename[recNum];
    //     const int type[3] = {0, 1, 2};
    //
    //     LOG << ENDL;
    //     for (int i = 0; i < 3; i++) {
    //         rangeCountFixWithLog<point>(wp, pkd, kdknn,
    //         singleQueryLogRepeatNum, type[i], recNum, Dim);
    //     }
    //     // rangeCountFix<point>(wp, pkd, kdknn, rounds, 2,
    //     rangeQueryNumInbaRatio, Dim);
    //     // rangeCountFix<point>(wp, pkd, kdknn, rounds, 2, recNum, Dim);
    //
    //     delete[] kdknn;
    // }
    //
    // if (queryType & (1 << 3)) {  // NOTE: range query
    //     if (summary == 0) {
    //         int recNum = rangeQueryNum;
    //         const int type[3] = {0, 1, 2};
    //
    //         LOG << ENDL;
    //         for (int i = 0; i < 3; i++) {
    //             //* run range count to obtain size
    //             kdknn = new Typename[recNum];
    //             points Out;
    //             rangeQuerySerialWithLog<point>(wp, pkd, kdknn,
    //             singleQueryLogRepeatNum, Out, type[i], recNum, Dim); delete[]
    //             kdknn;
    //         }
    //     } else if (summary == 1) {  // NOTE: for summary
    //         kdknn = new Typename[summaryRangeQueryNum];
    //         points Out;
    //         rangeQueryFix<point>(wp, pkd, kdknn, rounds, Out, 2,
    //         summaryRangeQueryNum, Dim); delete[] kdknn;
    //     }
    // }
    //
    // if (queryType & (1 << 4)) {  // NOTE: batch insertion with fraction
    //     const parlay::sequence<double> ratios = {0.0001, 0.0002, 0.0005,
    //     0.001, 0.002, 0.005, 0.01,
    //                                              0.02,   0.05,   0.1,    0.2,
    //                                              0.5,   1.0};
    //     for (int i = 0; i < ratios.size(); i++) {
    //         batchInsert<point>(pkd, wp, wi, Dim, rounds, ratios[i]);
    //     }
    // }
    //
    // if (queryType & (1 << 5)) {  // NOTE: batch deletion with fraction
    //     const parlay::sequence<double> ratios = {0.0001, 0.0002, 0.0005,
    //     0.001, 0.002, 0.005, 0.01,
    //                                              0.02,   0.05,   0.1,    0.2,
    //                                              0.5,   1.0};
    //     points tmp;
    //     for (int i = 0; i < ratios.size(); i++) {
    //         batchDelete<point>(pkd, wp, tmp, Dim, rounds, 0, ratios[i]);
    //     }
    // }
    //
    // if (queryType & (1 << 6)) {  //* incremental Build
    //     double step[4] = {0.1, 0.2, 0.25, 0.5};
    //     for (int i = 0; i < 4; i++) {
    //         incrementalBuild<point>(Dim, wp, rounds, pkd, step[i]);
    //     }
    // }
    //
    // if (queryType & (1 << 7)) {  //* incremental Delete
    //     double step[4] = {0.1, 0.2, 0.25, 0.5};
    //     for (int i = 0; i < 4; i++) {
    //         incrementalDelete<point>(Dim, wp, wi, rounds, pkd, step[i]);
    //     }
    // }
    //
    // if (queryType & (1 << 8)) {  //* batch insertion then knn
    //     kdknn = new Typename[wp.size()];
    //
    //     //* first normal build
    //     buildTree<point, 0>(Dim, wp, rounds, pkd);
    //     queryKNN<point>(Dim, wp, rounds, pkd, kdknn, K, false);
    //
    //     //* then incremental build
    //     incrementalBuild<point, 0>(Dim, wp, rounds, pkd, 0.1);
    //     queryKNN<point>(Dim, wp, rounds, pkd, kdknn, K, false);
    //
    //     delete[] kdknn;
    // }
    //
    // if (queryType & (1 << 9)) {  //* batch deletion then knn
    //     kdknn = new Typename[wp.size()];
    //
    //     //* first normal build
    //     buildTree<point, 0>(Dim, wp, rounds, pkd);
    //     queryKNN<point>(Dim, wp, rounds, pkd, kdknn, K, false);
    //
    //     //* then incremental delete
    //     incrementalDelete<point, 0>(Dim, wp, wi, rounds, pkd, 0.1);
    //     queryKNN<point>(Dim, wp, rounds, pkd, kdknn, K, false);
    //
    //     delete[] kdknn;
    // }
    //
    // if (queryType & (1 << 10)) {  // NOTE: test inbalance ratio
    //     const int fileNum = 10;
    //
    //     const size_t batchPointNum = wp.size() / fileNum;
    //
    //     points np, nq, up;
    //     std::string prefix, path;
    //     const string insertFileBack = insertFile;
    //     const string ten_varden_path =
    //     "/data/zmen002/kdtree/ss_varden/1000000000_3/10V.in"; const string
    //     one_uniform_nine_varden =
    //     "/data/zmen002/kdtree/ss_varden/1000000000_3/1U9V.in"; const string
    //     uniform_path = "/data/zmen002/kdtree/uniform/1000000000_3/1.in";
    //
    //     auto inbaQueryType = std::stoi(std::getenv("INBA_QUERY"));
    //     auto inbaBuildType = std::stoi(std::getenv("INBA_BUILD"));
    //
    //     // NOTE: helper functions
    //     auto clean = [&]() {
    //         prefix = insertFile.substr(0, insertFile.rfind("/"));
    //         np.clear();
    //         nq.clear();
    //     };
    //
    //     auto writeToFile = [&](string path) {
    //         std::ofstream f(path);
    //         f << np.size() << " " << Dim << std::endl;
    //         for (size_t i = 0; i < np.size(); i++) {
    //             for (size_t j = 0; j < Dim; j++) {
    //                 f << np[i].pnt[j] << " ";
    //             }
    //             f << std::endl;
    //         }
    //         f.close();
    //     };
    //
    //     // NOTE: run the test
    //     auto run = [&]() {
    //         if (inbaBuildType == 0) {
    //             buildTree<point, 2>(Dim, np, rounds, pkd);
    //         } else {
    //             // incrementalBuild<point, 2>(Dim, np, rounds, pkd,
    //             insertBatchInbaRatio);
    //
    //             size_t batchSize = static_cast<size_t>(up.size() *
    //             knnBatchInbaRatio); points newPts(batchSize);
    //             parlay::copy(up.cut(0, batchSize), newPts.cut(0, batchSize));
    //             incrementalBuildAndQuery<point, 2>(Dim, np, rounds, pkd,
    //             insertBatchInbaRatio, newPts);
    //         }
    //
    //         // if (inbaQueryType == 0) {
    //         //     size_t batchSize = static_cast<size_t>(np.size() *
    //         knnBatchInbaRatio);
    //         //     points newPts(batchSize);
    //         //     parlay::copy(np.cut(0, batchSize), newPts.cut(0,
    //         batchSize));
    //         //     kdknn = new Typename[batchSize];
    //         //     const int k[3] = {1, 5, 100};
    //         //     for (int i = 0; i < 3; i++) {
    //         //         queryKNN<point, 0, 1>(Dim, newPts, rounds, pkd, kdknn,
    //         k[i], true);
    //         //     }
    //         //     delete[] kdknn;
    //         // } else if (inbaQueryType == 1) {
    //         //     kdknn = new Typename[rangeQueryNumInbaRatio];
    //         //     int type = 2;
    //         //     rangeCountFix<point>(np, pkd, kdknn, rounds, type,
    //         rangeQueryNumInbaRatio, Dim);
    //         //     delete[] kdknn;
    //         // }
    //     };
    //
    //     read_points(uniform_path.c_str(), up, K);
    //
    //     LOG << "alpha: " << pkd.get_imbalance_ratio() << ENDL;
    //     // HACK: need start with varden file
    //     // NOTE: 1: 10*0.1 different vardens.
    //     clean();
    //     // for (int i = 1; i <= fileNum; i++) {
    //     //     path = prefix + "/" + std::to_string(i) + ".in";
    //     //     // std::cout << path << std::endl;
    //     //     read_points<point>(path.c_str(), nq, K);
    //     //     np.append(nq.cut(0, batchPointNum));
    //     //     nq.clear();
    //     // }
    //     // writeToFile(ten_varden_path);
    //     read_points(ten_varden_path.c_str(), np, K);
    //     assert(np.size() == wp.size());
    //     run();
    //
    //     // NOTE: 2: 1 uniform, and 9*0.1 same varden
    //     //* read varden first
    //     clean();
    //     // path = prefix + "/1.in";
    //     // std::cout << "varden path" << path << std::endl;
    //     // read_points<point>(path.c_str(), np, K);
    //     //* then read uniforprefixm
    //     // prefix = prefix.substr(0, prefix.rfind("/"));  // 1000000_3
    //     // prefix = prefix.substr(0, prefix.rfind("/"));  // ss_varden
    //     // path = prefix + "/uniform/" + std::to_string(wp.size()) + "_" +
    //     std::to_string(Dim) + "/1.in";
    //     // std::cout << "uniform path:" << path << std::endl;
    //
    //     // read_points<point>(path.c_str(), nq, K);
    //     // parlay::parallel_for(0, batchPointNum, [&](size_t i) { np[i] =
    //     nq[i]; });
    //     // writeToFile(one_uniform_nine_varden);
    //     read_points(one_uniform_nine_varden.c_str(), np, K);
    //     run();
    //
    //     //@ 3: 1 varden, but flatten;
    //     // clean();
    //     // path = prefix + "/1.in";
    //     // // std::cout << path << std::endl;
    //     // read_points<point>(path.c_str(), np, K);
    //     // buildTree<point, 0>(Dim, np, rounds, pkd);
    //     // pkd.flatten(pkd.get_root(), parlay::make_slice(np));
    //     // run();
    //
    //     // delete[] kdknn;
    // }
    //
    // if (queryType & (1 << 11)) {  // NOTE: osm by year
    //     // WARN: remember using double
    //     string osm_prefix = "/data/zmen002/kdtree/real_world/osm/year/";
    //     const std::vector<std::string> files = {"2014", "2015", "2016",
    //     "2017", "2018",
    //                                             "2019", "2020", "2021",
    //                                             "2022", "2023"};
    //     parlay::sequence<points> node_by_year(files.size());
    //     for (int i = 0; i < files.size(); i++) {
    //         std::string path = osm_prefix + "osm_" + files[i] + ".csv";
    //         // LOG << path << ENDL;
    //         read_points(path.c_str(), node_by_year[i], K);
    //     }
    //     kdknn = new Typename[batchQueryOsmSize];
    //     insertOsmByTime<point>(Dim, node_by_year, rounds, pkd, K, kdknn);
    //     delete[] kdknn;
    //
    //     // auto all_points = parlay::flatten(node_by_year);
    //     // queryKNN<point>(Dim, all_points, rounds, pkd, kdknn, K, false);
    // }
    //
    // if (queryType & (1 << 12)) {  // NOTE: osm by month
    //     // WARN: remember using double
    //     string osm_prefix = "/data/zmen002/kdtree/real_world/osm/month/";
    //     const std::vector<std::string> files = {"2014", "2015", "2016",
    //     "2017", "2018",
    //                                             "2019", "2020", "2021",
    //                                             "2022", "2023"};
    //     const std::vector<std::string> month = {"1", "2", "3", "4", "5", "6",
    //     "7", "8", "9", "10", "11", "12"};
    //
    //     parlay::sequence<points> node(files.size() * month.size());
    //     for (int i = 0; i < files.size(); i++) {
    //         for (int j = 0; j < month.size(); j++) {
    //             std::string path = osm_prefix + files[i] + "/" + month[j] +
    //             ".csv"; read_points(path.c_str(), node[i * month.size() + j],
    //             K);
    //         }
    //     }
    //     kdknn = new Typename[batchQueryOsmSize];
    //     insertOsmByTime<point>(Dim, node, rounds, pkd, K, kdknn);
    //     delete[] kdknn;
    //     // auto all_points = parlay::flatten(node);
    //     // queryKNN<point>(Dim, all_points, rounds, pkd, kdknn, K, false);
    // }
    //
    // if (queryType & (1 << 13)) {  // NOTE: serial insert VS batch insert
    //     // NOTE: first insert in serial one bu one
    //     const parlay::sequence<double> ratios = {1e-9, 2e-9, 5e-9, 1e-8,
    //     2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6,
    //                                              5e-6, 1e-5, 2e-5, 5e-5,
    //                                              1e-4, 2e-4, 5e-4, 1e-3,
    //                                              2e-3, 5e-3, 1e-2};
    //     // LOG << ENDL << "serial ";
    //     // batchInsert<point, true>(pkd, wp, wi, Dim, rounds,
    //     *ratios.rbegin()); LOG << ENDL; for (int i = 0; i < ratios.size();
    //     i++) {
    //         LOG << wi.size() * ratios[i] << " ";
    //         batchUpdateByStep<point, true>(pkd, wp, wi, Dim, rounds,
    //         ratios[i], *ratios.rbegin()); LOG << ENDL;
    //     }
    // }
    //
    // if (queryType & (1 << 14)) {  // NOTE: serial delete VS batch delete
    //     // NOTE: first insert in serial one bu one
    //     const parlay::sequence<double> ratios = {1e-9, 2e-9, 5e-9, 1e-8,
    //     2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6,
    //                                              5e-6, 1e-5, 2e-5, 5e-5,
    //                                              1e-4, 2e-4, 5e-4, 1e-3,
    //                                              2e-3, 5e-3, 1e-2};
    //     // LOG << ENDL << "serial ";
    //     // batchDelete<point, true>(pkd, wp, wi, Dim, rounds, false,
    //     *ratios.rbegin());
    //     // LOG << ENDL;
    //     for (int i = 0; i < ratios.size(); i++) {
    //         LOG << wi.size() * ratios[i] << " ";
    //         batchUpdateByStep<point, false>(pkd, wp, wp, Dim, rounds,
    //         ratios[i], *ratios.rbegin()); LOG << ENDL;
    //     }
    // }

    std::cout << std::endl << std::flush;

    tree.DeleteTree();

    return;
}

struct wrapper {
    struct QadTree {
        template<class Point>
        struct Desc {
            using TreeType = cpdd::QuadTree<Point, cpdd::RotateDim<Point>>;
        };
    };
    struct KDtree {
        template<class Point>
        struct Desc {
            using TreeType = cpdd::KdTree<Point, cpdd::MaxStretchDim<Point>>;
        };
    };
};

int main(int argc, char* argv[]) {
    commandLine P(argc, argv,
                  "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                  "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                  "<_insertFile>] [-s <summary>]");

    char* input_file_path = P.getOptionValue("-p");
    int K = P.getOptionIntValue("-k", 100);
    int Dim = P.getOptionIntValue("-d", 3);
    size_t N = P.getOptionLongValue("-n", -1);
    int tag = P.getOptionIntValue("-t", 1);
    int rounds = P.getOptionIntValue("-r", 3);
    int queryType = P.getOptionIntValue("-q", 0);
    int read_insert_file = P.getOptionIntValue("-i", 1);
    int summary = P.getOptionIntValue("-s", 0);
    int tree_type = P.getOptionIntValue("-T", 0);

    auto run_test = [&]<class Wrapper>(Wrapper) {
        auto run = [&](auto dim_wrapper) {
            constexpr const auto kDim = decltype(dim_wrapper)::value;
            using PointTypeAlias = PointType<Coord, kDim>;
            using Points = parlay::sequence<PointTypeAlias>;
            using Desc = typename Wrapper::template Desc<PointTypeAlias>;

            std::string name, insert_file_path = "";
            Points wp, wi;

            if (input_file_path != NULL) {  // NOTE: read main points
                name = std::string(input_file_path);
                name = name.substr(name.rfind("/") + 1);
                std::cout << name << " ";
                auto [n, d] =
                    read_points<PointTypeAlias>(input_file_path, wp, K);
                N = n;
                assert(d == Dim);
            }

            if (read_insert_file == 1) {  // NOTE: read points to be inserted
                int id = std::stoi(name.substr(0, name.find_first_of('.')));
                id = (id + 1) % 3;  // WARN: MOD graph number used to test
                if (!id) id++;
                int pos = std::string(input_file_path).rfind("/") + 1;
                insert_file_path = std::string(input_file_path).substr(0, pos) +
                                   std::to_string(id) + ".in";
                auto [n, d] = read_points<PointTypeAlias>(
                    insert_file_path.c_str(), wi, K);
                assert(d == Dim);
            }

            TestSpacialTree<Desc>(Dim, wp, wi, N, K, rounds, insert_file_path,
                                  tag, queryType, summary);
        };

        if (tag == -1) {
            // NOTE: serial run
            ;
        } else if (Dim == 2) {
            run(std::integral_constant<int, 2>{});
        } else if (Dim == 3) {
            run(std::integral_constant<int, 3>{});
        } else if (Dim == 5) {
            run(std::integral_constant<int, 5>{});
        } else if (Dim == 7) {
            run(std::integral_constant<int, 7>{});
        } else if (Dim == 9) {
            run(std::integral_constant<int, 9>{});
        } else if (Dim == 10) {
            run(std::integral_constant<int, 10>{});
        }
    };

    if (tree_type == 0) {
        run_test(wrapper::KDtree{});
    } else if (tree_type == 1) {
        run_test(wrapper::QadTree{});
    } else if (tree_type == 2) {
    }

    return 0;
}
