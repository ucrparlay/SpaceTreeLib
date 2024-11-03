
#include <algorithm>
#include <cstddef>
#include <cstdlib>

#include "cpdd/dependence/splitter.h"
#include "cpdd/kd_tree.h"
#include "parlay/primitives.h"
#include "testFramework.h"

template <class TreeDesc, typename Point>
void TestSpacialTree(int const& kDim, parlay::sequence<Point> const& wp,
                     parlay::sequence<Point> const& wi,
                     [[maybe_unused]] size_t const& N, int const& K,
                     int const& kRounds,
                     [[maybe_unused]] string const& kInsertFile,
                     int const& kTag, int const& kQueryType,
                     int const kSummary) {
  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;

  Tree tree;

  buildTree<Point, Tree, 2>(kDim, wp, kRounds, tree);

  Typename* kdknn = nullptr;

  // NOTE: batch insert
  if (kTag & (1 << 0)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchInsert<Point, Tree, true>(tree, wp, wi, kDim, kRounds, ratios[i]);
      }
    } else {
      BatchInsert<Point, Tree, true>(tree, wp, wi, kDim, kRounds,
                                     kBatchInsertRatio);
    }
  }

  // NOTE: batch delete
  if (kTag & (1 << 1)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchDelete<Point, Tree, true>(tree, wp, wi, kDim, kRounds, 0,
                                       ratios[i]);
      }
    } else {
      BatchDelete<Point, Tree, true>(tree, wp, wi, kDim, kRounds, 0,
                                     kBatchInsertRatio);
    }
  }

  // NOTE: batch diff
  if (kTag & (1 << 2)) {
    if (kSummary) {
      parlay::sequence<double> const overlap_ratios = {0.1, 0.2, 0.5, 1};
      for (size_t i = 0; i < overlap_ratios.size(); i++) {
        BatchDiff<Point, Tree, true>(tree, wp, kDim, kRounds,
                                     kBatchDiffTotalRatio, overlap_ratios[i]);
      }
    } else {
      BatchDiff<Point, Tree, true>(tree, wp, kDim, kRounds,
                                   kBatchDiffTotalRatio,
                                   kBatchDiffOverlapRatio);
    }
  }

  if (kQueryType & (1 << 0)) {  // NOTE: KNN
    auto run_batch_knn = [&](Points const& pts, int kth, size_t batchSize) {
      Points newPts(batchSize);
      parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
      kdknn = new Typename[batchSize];
      queryKNN<Point>(kDim, newPts, kRounds, tree, kdknn, kth, true);
      delete[] kdknn;
    };

    size_t batchSize = static_cast<size_t>(wp.size() * batchQueryRatio);

    if (kSummary == 0) {
      int k[3] = {1, 10, 100};
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp, k[i], batchSize);
      }
    } else {  // test kSummary
      run_batch_knn(wp, K, batchSize);
    }
  }

  if (kQueryType & (1 << 1)) {  // NOTE: batch NN query

    auto run_batch_knn = [&](Points const& pts, size_t batchSize) {
      Points newPts(batchSize);
      parlay::copy(pts.cut(0, batchSize), newPts.cut(0, batchSize));
      kdknn = new Typename[batchSize];
      queryKNN<Point, Tree, true, true>(kDim, newPts, kRounds, tree, kdknn, K,
                                        true);
      delete[] kdknn;
    };

    run_batch_knn(wp, static_cast<size_t>(wp.size() * batchQueryRatio));
    std::vector<double> const batchRatios = {0.001, 0.01, 0.1, 0.2, 0.5};
    for (auto ratio : batchRatios) {
      run_batch_knn(wp, static_cast<size_t>(wp.size() * ratio));
    }
    // for (auto ratio : batchRatios) {
    //   run_batch_knn(wi, static_cast<size_t>(wi.size() * ratio));
    // }
  }

  if (kQueryType & (1 << 2)) {  // NOTE: range count
    int recNum = rangeQueryNum;
    kdknn = new Typename[recNum];
    // int const type[3] = {0, 1, 2};

    // std::cout << std::endl;
    for (int i = 0; i < 3; i++) {
      rangeCountFix<Point>(wp, tree, kdknn, kRounds, i, recNum, kDim);
      // rangeCountFixWithLog<Point>(wp, tree, kdknn,
      //                             singleQueryLogRepeatNum, type[i],
      //                             recNum, kDim);
    }
    // rangeCountFix<Point>(wp, tree, kdknn, kRounds, 2, recNum, kDim);

    delete[] kdknn;
  }

  if (kQueryType & (1 << 3)) {  // NOTE: range query
    if (kSummary == 0) {
      int recNum = rangeQueryNum;
      kdknn = new Typename[recNum];
      // int const type[3] = {0, 1, 2};

      // std::cout << std::endl;
      for (int i = 0; i < 3; i++) {
        //* run range count to obtain size
        Points Out;
        rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, i, recNum, kDim);
        // rangeQuerySerialWithLog<Point>(wp, tree, kdknn,
        //                                singleQueryLogRepeatNum, Out,
        //                                type[i], recNum, kDim);
      }
      delete[] kdknn;
    } else if (kSummary == 1) {  // NOTE: for kSummary
      kdknn = new Typename[summaryRangeQueryNum];
      Points Out;
      rangeQueryFix<Point>(wp, tree, kdknn, kRounds, Out, 2,
                           summaryRangeQueryNum, kDim);
      delete[] kdknn;
    }
  }
  //
  // if (kQueryType & (1 << 4)) {  // NOTE: batch insertion with fraction
  //     const parlay::sequence<double> ratios = {0.0001, 0.0002, 0.0005,
  //     0.001, 0.002, 0.005, 0.01,
  //                                              0.02,   0.05,   0.1,    0.2,
  //                                              0.5,   1.0};
  //     for (int i = 0; i < ratios.size(); i++) {
  //         BatchInsert<Point>(tree, wp, wi, kDim, kRounds, ratios[i]);
  //     }
  // }
  //
  // if (kQueryType & (1 << 5)) {  // NOTE: batch deletion with fraction
  //     const parlay::sequence<double> ratios = {0.0001, 0.0002, 0.0005,
  //     0.001, 0.002, 0.005, 0.01,
  //                                              0.02,   0.05,   0.1,    0.2,
  //                                              0.5,   1.0};
  //     Points tmp;
  //     for (int i = 0; i < ratios.size(); i++) {
  //         batchDelete<Point>(tree, wp, tmp, kDim, kRounds, 0, ratios[i]);
  //     }
  // }
  //
  // if (kQueryType & (1 << 6)) {  //* incremental Build
  //     double step[4] = {0.1, 0.2, 0.25, 0.5};
  //     for (int i = 0; i < 4; i++) {
  //         incrementalBuild<Point>(kDim, wp, kRounds, tree, step[i]);
  //     }
  // }
  //
  // if (kQueryType & (1 << 7)) {  //* incremental Delete
  //     double step[4] = {0.1, 0.2, 0.25, 0.5};
  //     for (int i = 0; i < 4; i++) {
  //         incrementalDelete<Point>(kDim, wp, wi, kRounds, tree, step[i]);
  //     }
  // }
  //
  // if (kQueryType & (1 << 8)) {  //* batch insertion then knn
  //     kdknn = new Typename[wp.size()];
  //
  //     //* first normal build
  //     buildTree<Point, 0>(kDim, wp, kRounds, tree);
  //     queryKNN<Point>(kDim, wp, kRounds, tree, kdknn, K, false);
  //
  //     //* then incremental build
  //     incrementalBuild<Point, 0>(kDim, wp, kRounds, tree, 0.1);
  //     queryKNN<Point>(kDim, wp, kRounds, tree, kdknn, K, false);
  //
  //     delete[] kdknn;
  // }
  //
  // if (kQueryType & (1 << 9)) {  //* batch deletion then knn
  //     kdknn = new Typename[wp.size()];
  //
  //     //* first normal build
  //     buildTree<Point, 0>(kDim, wp, kRounds, tree);
  //     queryKNN<Point>(kDim, wp, kRounds, tree, kdknn, K, false);
  //
  //     //* then incremental delete
  //     incrementalDelete<Point, 0>(kDim, wp, wi, kRounds, tree, 0.1);
  //     queryKNN<Point>(kDim, wp, kRounds, tree, kdknn, K, false);
  //
  //     delete[] kdknn;
  // }
  //
  // if (kQueryType & (1 << 10)) {  // NOTE: test inbalance ratio
  //     const int fileNum = 10;
  //
  //     const size_t batchPointNum = wp.size() / fileNum;
  //
  //     Points np, nq, up;
  //     std::string prefix, path;
  //     const string insertFileBack = insertFile;
  //     const string ten_varden_path =
  //     "/data/zmen002/kdtree/ss_varden/1000000000_3/10V.in"; const string
  //     one_uniform_nine_varden =
  //     "/data/zmen002/kdtree/ss_varden/1000000000_3/1U9V.in"; const string
  //     uniform_path = "/data/zmen002/kdtree/uniform/1000000000_3/1.in";
  //
  //     auto inbakQueryType = std::stoi(std::getenv("INBA_QUERY"));
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
  //         f << np.size() << " " << kDim << std::endl;
  //         for (size_t i = 0; i < np.size(); i++) {
  //             for (size_t j = 0; j < kDim; j++) {
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
  //             buildTree<Point, 2>(kDim, np, kRounds, tree);
  //         } else {
  //             // incrementalBuild<Point, 2>(kDim, np, kRounds, tree,
  //             insertBatchInbaRatio);
  //
  //             size_t batchSize = static_cast<size_t>(up.size() *
  //             knnBatchInbaRatio); Points newPts(batchSize);
  //             parlay::copy(up.cut(0, batchSize), newPts.cut(0, batchSize));
  //             incrementalBuildAndQuery<Point, 2>(kDim, np, kRounds, tree,
  //             insertBatchInbaRatio, newPts);
  //         }
  //
  //         // if (inbakQueryType == 0) {
  //         //     size_t batchSize = static_cast<size_t>(np.size() *
  //         knnBatchInbaRatio);
  //         //     Points newPts(batchSize);
  //         //     parlay::copy(np.cut(0, batchSize), newPts.cut(0,
  //         batchSize));
  //         //     kdknn = new Typename[batchSize];
  //         //     const int k[3] = {1, 5, 100};
  //         //     for (int i = 0; i < 3; i++) {
  //         //         queryKNN<Point, 0, 1>(kDim, newPts, kRounds, tree,
  //         kdknn, k[i], true);
  //         //     }
  //         //     delete[] kdknn;
  //         // } else if (inbakQueryType == 1) {
  //         //     kdknn = new Typename[rangeQueryNumInbaRatio];
  //         //     int type = 2;
  //         //     rangeCountFix<Point>(np, tree, kdknn, kRounds, type,
  //         rangeQueryNumInbaRatio, kDim);
  //         //     delete[] kdknn;
  //         // }
  //     };
  //
  //     read_points(uniform_path.c_str(), up, K);
  //
  //     std::cout << "alpha: " << tree.get_imbalance_ratio() << std::endl;
  //     // HACK: need start with varden file
  //     // NOTE: 1: 10*0.1 different vardens.
  //     clean();
  //     // for (int i = 1; i <= fileNum; i++) {
  //     //     path = prefix + "/" + std::to_string(i) + ".in";
  //     //     // std::cout << path << std::endl;
  //     //     read_points<Point>(path.c_str(), nq, K);
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
  //     // read_points<Point>(path.c_str(), np, K);
  //     //* then read uniforprefixm
  //     // prefix = prefix.substr(0, prefix.rfind("/"));  // 1000000_3
  //     // prefix = prefix.substr(0, prefix.rfind("/"));  // ss_varden
  //     // path = prefix + "/uniform/" + std::to_string(wp.size()) + "_" +
  //     std::to_string(kDim) + "/1.in";
  //     // std::cout << "uniform path:" << path << std::endl;
  //
  //     // read_points<Point>(path.c_str(), nq, K);
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
  //     // read_points<Point>(path.c_str(), np, K);
  //     // buildTree<Point, 0>(kDim, np, kRounds, tree);
  //     // tree.flatten(tree.get_root(), parlay::make_slice(np));
  //     // run();
  //
  //     // delete[] kdknn;
  // }
  //
  // if (kQueryType & (1 << 11)) {  // NOTE: osm by year
  //     // WARN: remember using double
  //     string osm_prefix = "/data/zmen002/kdtree/real_world/osm/year/";
  //     const std::vector<std::string> files = {"2014", "2015", "2016",
  //     "2017", "2018",
  //                                             "2019", "2020", "2021",
  //                                             "2022", "2023"};
  //     parlay::sequence<Points> node_by_year(files.size());
  //     for (int i = 0; i < files.size(); i++) {
  //         std::string path = osm_prefix + "osm_" + files[i] + ".csv";
  //         // std::cout << path << std::endl;
  //         read_points(path.c_str(), node_by_year[i], K);
  //     }
  //     kdknn = new Typename[batchQueryOsmSize];
  //     insertOsmByTime<Point>(kDim, node_by_year, kRounds, tree, K, kdknn);
  //     delete[] kdknn;
  //
  //     // auto all_points = parlay::flatten(node_by_year);
  //     // queryKNN<Point>(kDim, all_points, kRounds, tree, kdknn, K, false);
  // }
  //
  // if (kQueryType & (1 << 12)) {  // NOTE: osm by month
  //     // WARN: remember using double
  //     string osm_prefix = "/data/zmen002/kdtree/real_world/osm/month/";
  //     const std::vector<std::string> files = {"2014", "2015", "2016",
  //     "2017", "2018",
  //                                             "2019", "2020", "2021",
  //                                             "2022", "2023"};
  //     const std::vector<std::string> month = {"1", "2", "3", "4", "5", "6",
  //     "7", "8", "9", "10", "11", "12"};
  //
  //     parlay::sequence<Points> node(files.size() * month.size());
  //     for (int i = 0; i < files.size(); i++) {
  //         for (int j = 0; j < month.size(); j++) {
  //             std::string path = osm_prefix + files[i] + "/" + month[j] +
  //             ".csv"; read_points(path.c_str(), node[i * month.size() + j],
  //             K);
  //         }
  //     }
  //     kdknn = new Typename[batchQueryOsmSize];
  //     insertOsmByTime<Point>(kDim, node, kRounds, tree, K, kdknn);
  //     delete[] kdknn;
  //     // auto all_points = parlay::flatten(node);
  //     // queryKNN<Point>(kDim, all_points, kRounds, tree, kdknn, K, false);
  // }
  //
  // if (kQueryType & (1 << 13)) {  // NOTE: serial insert VS batch insert
  //     // NOTE: first insert in serial one bu one
  //     const parlay::sequence<double> ratios = {1e-9, 2e-9, 5e-9, 1e-8,
  //     2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6,
  //                                              5e-6, 1e-5, 2e-5, 5e-5,
  //                                              1e-4, 2e-4, 5e-4, 1e-3,
  //                                              2e-3, 5e-3, 1e-2};
  //     // std::cout << std::endl << "serial ";
  //     // BatchInsert<Point, true>(tree, wp, wi, kDim, kRounds,
  //     *ratios.rbegin()); std::cout << std::endl; for (int i = 0; i <
  //     ratios.size(); i++) {
  //         std::cout << wi.size() * ratios[i] << " ";
  //         batchUpdateByStep<Point, true>(tree, wp, wi, kDim, kRounds,
  //         ratios[i], *ratios.rbegin()); std::cout << std::endl;
  //     }
  // }
  //
  // if (kQueryType & (1 << 14)) {  // NOTE: serial delete VS batch delete
  //     // NOTE: first insert in serial one bu one
  //     const parlay::sequence<double> ratios = {1e-9, 2e-9, 5e-9, 1e-8,
  //     2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6,
  //                                              5e-6, 1e-5, 2e-5, 5e-5,
  //                                              1e-4, 2e-4, 5e-4, 1e-3,
  //                                              2e-3, 5e-3, 1e-2};
  //     // std::cout << std::endl << "serial ";
  //     // batchDelete<Point, true>(tree, wp, wi, kDim, kRounds, false,
  //     *ratios.rbegin());
  //     // std::cout << std::endl;
  //     for (int i = 0; i < ratios.size(); i++) {
  //         std::cout << wi.size() * ratios[i] << " ";
  //         batchUpdateByStep<Point, false>(tree, wp, wp, kDim, kRounds,
  //         ratios[i], *ratios.rbegin()); std::cout << std::endl;
  //     }
  // }

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
  int tag = P.getOptionIntValue("-t", 1);
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
        auto [n, d] =
            read_points<PointTypeAlias>(insert_file_path.c_str(), wi, K);
        assert(d == kDim);
      }

      TestSpacialTree<Desc>(kDim, wp, wi, N, K, kRounds, insert_file_path, tag,
                            kQueryType, kSummary);
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
    run_test(wrapper::KDtree{});
  } else if (tree_type == 1 && kDim == 2) {
    // run_test(wrapper::QuadTree{});
  }
  // else if (tree_type == 1 && kDim == 3) {
  //     run_test(wrapper::OctTree{});
  // }

  return 0;
}
