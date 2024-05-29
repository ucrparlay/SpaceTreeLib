#include <algorithm>
#include <cstdlib>
#include "testFramework.h"

template<typename Point>
void testParallelKDtree(const int& Dim, const int& kLeaveWrap, parlay::sequence<Point>& wp, const size_t& N,
                        const int& K, const int& rounds, const string& insertFile, const int& tag,
                        const int& queryType) {
    using tree = BaseTree<Point>;
    using Points = typename tree::Points;
    using node = typename tree::node;
    using interior = typename tree::interior;
    using leaf = typename tree::leaf;
    using NodeTag = typename tree::NodeTag;
    using NodeTagSeq = typename tree::NodeTagSeq;
    using Box = typename tree::Box;
    using boxs = parlay::sequence<Box>;

    // auto boxs = gen_rectangles<Point>( 1000, 2, wp, Dim );
    // for ( int i = 0; i < 10; i++ ) {
    //   LOG << boxs[i].first << " " << boxs[i].second << ENDL;
    // }
    // return;

    if (N != wp.size()) {
        puts("input parameter N is different to input Points size");
        abort();
    }

    tree pkd;

    Points wi;
    if (insertFile != "") {
        auto [nn, nd] = read_points<Point>(insertFile.c_str(), wi, K);
        if (nd != Dim) {
            puts("read inserted Points dimension wrong");
            abort();
        }
    }

    Typename* kdknn = nullptr;

    //* begin test
    buildTree<Point, tree>(Dim, wp, rounds, pkd);

    //* batch insert
    if (tag >= 1) {
        BatchInsert<Point>(pkd, wp, wi, Dim, rounds, batchInsertRatio);
    }

    //* batch delete
    if (tag >= 2) {
        BatchDelete<Point>(pkd, wp, wi, Dim, rounds, 0, batchInsertRatio);
    }

    if (queryType & (1 << 0)) {  //* KNN
        kdknn = new Typename[wp.size()];
        if (tag == 0) {
            int k[3] = {1, 10, 100};
            for (int i = 0; i < 3; i++) {
                queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, k[i], false);
                // veryLargeKNN<Point>( Dim, wp, rounds, pkd, kdknn, k[i], false );
            }
        } else {  // test summary
            queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, K, false);
        }
        delete[] kdknn;
    }

    if (queryType & (1 << 1)) {  //* batch NN query
        Points new_wp(batchQuerySize);
        parlay::copy(wp.cut(0, batchQuerySize), new_wp.cut(0, batchQuerySize));
        kdknn = new Typename[batchQuerySize];

        queryKNN(Dim, new_wp, rounds, pkd, kdknn, K, true);

        delete[] kdknn;
    }

    int recNum = rangeQueryNum;

    if (queryType & (1 << 2)) {  //* range count
        kdknn = new Typename[recNum];
        int type[3] = {0, 1, 2};

        for (int i = 0; i < 3; i++) {
            rangeCountFix<Point>(wp, pkd, kdknn, rounds, type[i], recNum, Dim);
        }

        delete[] kdknn;
    }

    if (queryType & (1 << 3)) {  //* range query

        if (tag == 0) {
            const int type[3] = {0, 1, 2};
            for (int i = 0; i < 3; i++) {
                //* run range count to obtain size
                kdknn = new Typename[recNum];
                Points Out;
                //* range query
                rangeQueryFix<Point>(wp, pkd, kdknn, rounds, Out, type[i], recNum, Dim);
            }
        } else if (tag == 2) {
            kdknn = new Typename[recNum];
            Points Out;
            rangeQueryFix<Point>(wp, pkd, kdknn, rounds, Out, 2, recNum, Dim);
        }

        delete[] kdknn;
    }

    // if ( queryType & ( 1 << 3 ) ) {  //* generate knn
    // generate_knn<Point>( Dim, wp, K,
    // "/data9/zmen002/knn/GeoLifeNoScale.pbbs.out"
    // );
    // }

    if (queryType & (1 << 4)) {  //* batch insertion with fraction
        double ratios[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        for (int i = 0; i < 10; i++) {
            BatchInsert<Point>(pkd, wp, wi, Dim, rounds, ratios[i]);
        }
    }

    if (queryType & (1 << 5)) {  //* batch deletion with fraction
        const double ratios[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        // double ratios[10] = { 1.0 };
        Points tmp;
        for (int i = 0; i < 10; i++) {
            BatchDelete<Point>(pkd, wp, tmp, Dim, rounds, 0, ratios[i]);
        }
    }

    if (queryType & (1 << 6)) {  //* incremental Build
        double step[4] = {0.1, 0.2, 0.25, 0.5};
        for (int i = 0; i < 4; i++) {
            incrementalBuild<Point>(Dim, wp, rounds, pkd, step[i]);
        }
    }

    if (queryType & (1 << 7)) {  //* incremental Delete
        double step[4] = {0.1, 0.2, 0.25, 0.5};
        for (int i = 0; i < 4; i++) {
            incrementalDelete<Point>(Dim, wp, wi, rounds, pkd, step[i]);
        }
    }

    if (queryType & (1 << 8)) {  //* batch insertion then knn
        kdknn = new Typename[wp.size()];

        //* first normal Build
        buildTree<Point, tree, 0>(Dim, wp, rounds, pkd);
        queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, K, false);

        //* then incremental Build
        incrementalBuild<Point, 0>(Dim, wp, rounds, pkd, 0.1);
        queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, K, false);

        delete[] kdknn;
    }

    if (queryType & (1 << 9)) {  //* batch deletion then knn
        kdknn = new Typename[wp.size()];

        //* first normal Build
        buildTree<Point, tree, 0>(Dim, wp, rounds, pkd);
        queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, K, false);

        //* then incremental delete
        incrementalDelete<Point, 0>(Dim, wp, wi, rounds, pkd, 0.1);
        queryKNN<Point>(Dim, wp, rounds, pkd, kdknn, K, false);

        delete[] kdknn;
    }

    if (queryType & (1 << 10)) {  //* test inbalance ratio
        Points np, nq;
        std::string prefix, path;
        size_t batchSize = wp.size() / 10;
        kdknn = new Typename[wp.size()];

        auto inbaQueryType = std::stoi(std::getenv("INBA_QUERY"));
        auto inbaBuildType = std::stoi(std::getenv("INBA_BUILD"));

        //@ helper functions
        auto clean = [&]() {
            prefix = insertFile.substr(0, insertFile.rfind("/"));
            np.clear();
            nq.clear();
        };

        auto run = [&]() {
            if (inbaBuildType == 0) {
                buildTree<Point, tree, 2>(Dim, np, rounds, pkd);
            } else {
                incrementalBuild<Point, 2>(Dim, np, rounds, pkd, 0.01);
            }

            if (inbaQueryType == 0) {
                const int k[3] = {1, 5, 100};
                for (int i = 0; i < 3; i++) {
                    queryKNN<Point, 0, 1>(Dim, np, rounds, pkd, kdknn, k[i], false);
                }
            } else if (inbaQueryType == 1) {
                int type = 2;
                rangeCountFix<Point>(wp, pkd, kdknn, rounds, type, rangeQueryNumInbaRatio, Dim);
            }
        };

        LOG << "alpha: " << pkd.GetImbalanceRatio() << ENDL;
        //! start with varden file
        //@ 1: 10*0.1 different vardens.
        clean();
        for (int i = 1; i <= 10; i++) {
            path = prefix + "/" + std::to_string(i) + ".in";
            // std::cout << path << std::endl;
            read_points<Point>(path.c_str(), nq, K);
            np.append(nq.cut(0, batchSize));
            nq.clear();
        }
        assert(np.size() == wp.size());
        run();

        //@ 2: 1 uniform, and 9*0.1 same varden
        //* read varden first
        clean();
        path = prefix + "/1.in";
        // std::cout << path << std::endl;
        read_points<Point>(path.c_str(), np, K);
        //* then read uniforprefixm
        prefix = prefix.substr(0, prefix.rfind("/"));  // 1000000_3
        prefix = prefix.substr(0, prefix.rfind("/"));  // ss_varden
        path = prefix + "/uniform/" + std::to_string(np.size()) + "_" + std::to_string(Dim) + "/1.in";
        // std::cout << path << std::endl;

        read_points<Point>(path.c_str(), nq, K);
        parlay::parallel_for(0, batchSize, [&](size_t i) { np[i] = nq[i]; });
        run();

        //@ 3: 1 varden, but Flatten;
        clean();
        path = prefix + "/1.in";
        // std::cout << path << std::endl;
        read_points<Point>(path.c_str(), np, K);
        buildTree<Point, tree, 0>(Dim, np, rounds, pkd);
        pkd.Flatten(pkd.GetRoot(), parlay::make_slice(np));
        run();

        delete[] kdknn;
    }

    // generate_knn( Dim, wp, K, "knn.out" );

    std::cout << std::endl << std::flush;

    return;
}

int main(int argc, char* argv[]) {
    commandLine P(argc, argv,
                  "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                  "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                  "<_insertFile>]");
    char* iFile = P.getOptionValue("-p");
    int K = P.getOptionIntValue("-k", 100);
    int Dim = P.getOptionIntValue("-d", 3);
    size_t N = P.getOptionLongValue("-n", -1);
    int tag = P.getOptionIntValue("-t", 1);
    int rounds = P.getOptionIntValue("-r", 3);
    int queryType = P.getOptionIntValue("-q", 0);
    int readInsertFile = P.getOptionIntValue("-i", 1);

    int kLeaveWrap = 32;
    parlay::sequence<PointType<Coord, 15>> wp;
    // parlay::sequence<PointID<Coord, 15>> wp;
    std::string name, insertFile = "";

    //* initialize Points
    if (iFile != NULL) {
        name = std::string(iFile);
        name = name.substr(name.rfind("/") + 1);
        std::cout << name << " ";
        auto [n, d] = read_points<PointType<Coord, 15>>(iFile, wp, K);
        // auto [n, d] = read_points<PointID<Coord, 15>>( iFile, wp, K );
        N = n;
        assert(d == Dim);
    } else {  //* construct data byself
        K = 100;
        generate_random_points<PointType<Coord, 15>>(wp, 1000000, N, Dim);
        // generate_random_points<PointID<Coord, 15>>( wp, 1000000, N, Dim );
        std::string name = std::to_string(N) + "_" + std::to_string(Dim) + ".in";
        std::cout << name << " ";
    }

    if (readInsertFile == 1) {
        int id = std::stoi(name.substr(0, name.find_first_of('.')));
        id = (id + 1) % 3;  //! MOD graph number used to test
        if (!id) id++;
        int pos = std::string(iFile).rfind("/") + 1;
        insertFile = std::string(iFile).substr(0, pos) + std::to_string(id) + ".in";
    }

    assert(N > 0 && Dim > 0 && K > 0 && kLeaveWrap >= 1);

    // if ( tag == -1 ) {
    //   //* serial run
    //   // todo rewrite test serial code
    //   // testSerialKDtree( Dim, kLeaveWrap, wp, N, K );
    // } else if ( Dim == 2 ) {
    //   auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointID<Coord, 2> {
    //     return PointID<Coord, 2>( wp[i].pnt.begin(), i );
    //   } );
    //   decltype( wp )().swap( wp );
    //   testParallelKDtree<PointID<Coord, 2>>( Dim, kLeaveWrap, pts, N, K,
    //   rounds, insertFile,
    //                                          tag, queryType );
    // } else if ( Dim == 3 ) {
    //   auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointID<Coord, 3> {
    //     return PointID<Coord, 3>( wp[i].pnt.begin(), i );
    //   } );
    //   decltype( wp )().swap( wp );
    //   testParallelKDtree<PointID<Coord, 3>>( Dim, kLeaveWrap, pts, N, K,
    //   rounds, insertFile,
    //                                          tag, queryType );
    // } else if ( Dim == 5 ) {
    //   auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointID<Coord, 5> {
    //     return PointID<Coord, 5>( wp[i].pnt.begin(), i );
    //   } );
    //   decltype( wp )().swap( wp );
    //   testParallelKDtree<PointID<Coord, 5>>( Dim, kLeaveWrap, pts, N, K,
    //   rounds, insertFile,
    //                                          tag, queryType );
    // } else if ( Dim == 7 ) {
    //   auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointID<Coord, 7> {
    //     return PointID<Coord, 7>( wp[i].pnt.begin(), i );
    //   } );
    //   decltype( wp )().swap( wp );
    //   testParallelKDtree<PointID<Coord, 7>>( Dim, kLeaveWrap, pts, N, K,
    //   rounds, insertFile,
    //                                          tag, queryType );
    // }
    //
    if (tag == -1) {
        //* serial run
        // todo rewrite test serial code
        // testSerialKDtree( Dim, kLeaveWrap, wp, N, K );
    } else if (Dim == 2) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 2> { return PointType<Coord, 2>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 2>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 3) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 3> { return PointType<Coord, 3>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 3>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 5) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 5> { return PointType<Coord, 5>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 5>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 7) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 7> { return PointType<Coord, 7>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 7>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 9) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 9> { return PointType<Coord, 9>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 9>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 10) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 10> { return PointType<Coord, 10>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testParallelKDtree<PointType<Coord, 10>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    }

    return 0;
}
