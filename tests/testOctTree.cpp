#include <algorithm>
#include <cstdlib>
#include "testFramework.h"

template<typename Point>
void testOctTree(const int& Dim, const int& kLeaveWrap, parlay::sequence<Point>& wp, const size_t& N, const int& K,
                 const int& rounds, const string& insertFile, const int& tag, const int& queryType) {
    using tree = octTree<Point>;
    tree pkd;
    buildTree<Point, tree>(Dim, wp, rounds, pkd);

    std::cout << std::endl << std::flush;
    pkd.DeleteTree();
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
        // auto [n, d] = read_points<PointID<Coord, 15>>(iFile, wp, K);
        N = n;
        assert(d == Dim);
    } else {  //* construct data byself
        K = 100;
        generate_random_points<PointType<Coord, 15>>(wp, 1000000, N, Dim);
        // generate_random_points<PointID<Coord, 15>>(wp, 1000000, N, Dim);
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

    if (tag == -1) {
        //* serial run
        // todo rewrite test serial code
        // testSerialKDtree( Dim, kLeaveWrap, wp, N, K );
    } else if (Dim == 2) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 2> { return PointType<Coord, 2>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 2>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 3) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 3> { return PointType<Coord, 3>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 3>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 5) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 5> { return PointType<Coord, 5>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 5>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 7) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 7> { return PointType<Coord, 7>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 7>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 9) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 9> { return PointType<Coord, 9>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 9>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    } else if (Dim == 10) {
        auto pts = parlay::tabulate(
            N, [&](size_t i) -> PointType<Coord, 10> { return PointType<Coord, 10>(wp[i].pnt.begin()); });
        decltype(wp)().swap(wp);
        testOctTree<PointType<Coord, 10>>(Dim, kLeaveWrap, pts, N, K, rounds, insertFile, tag, queryType);
    }
    // if (tag == -1) {
    //   //* serial run
    //   // todo rewrite test serial code
    //   // testSerialKDtree( Dim, kLeaveWrap, wp, N, K );
    // } else if (Dim == 2) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 2> {
    //     return PointID<Coord, 2>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 2>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                  insertFile, tag, queryType);
    // } else if (Dim == 3) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 3> {
    //     return PointID<Coord, 3>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 3>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                  insertFile, tag, queryType);
    // } else if (Dim == 5) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 5> {
    //     return PointID<Coord, 5>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 5>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                  insertFile, tag, queryType);
    // } else if (Dim == 7) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 7> {
    //     return PointID<Coord, 7>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 7>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                  insertFile, tag, queryType);
    // } else if (Dim == 9) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 9> {
    //     return PointID<Coord, 9>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 9>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                  insertFile, tag, queryType);
    // } else if (Dim == 10) {
    //   auto pts = parlay::tabulate(N, [&](size_t i) -> PointID<Coord, 10> {
    //     return PointID<Coord, 10>(wp[i].pnt.begin(), wp[i].get_id());
    //   });
    //   decltype(wp)().swap(wp);
    //   testOctTree<PointID<Coord, 10>>(Dim, kLeaveWrap, pts, N, K, rounds,
    //                                   insertFile, tag, queryType);
    // }

    return 0;
}
