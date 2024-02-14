#include <algorithm>
#include <cstdlib>
#include "testFramework.h"

template<typename point>
void
testOctTree( const int& Dim, const int& LEAVE_WRAP, parlay::sequence<point>& wp,
             const size_t& N, const int& K, const int& rounds, const string& insertFile,
             const int& tag, const int& queryType ) {
  using tree = octTree<point>;
  tree pkd;
  buildTree<point, tree>( Dim, wp, rounds, pkd );

  std::cout << std::endl << std::flush;

  return;
}

int
main( int argc, char* argv[] ) {
  commandLine P( argc, argv,
                 "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                 "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                 "<_insertFile>]" );
  char* iFile = P.getOptionValue( "-p" );
  int K = P.getOptionIntValue( "-k", 100 );
  int Dim = P.getOptionIntValue( "-d", 3 );
  size_t N = P.getOptionLongValue( "-n", -1 );
  int tag = P.getOptionIntValue( "-t", 1 );
  int rounds = P.getOptionIntValue( "-r", 3 );
  int queryType = P.getOptionIntValue( "-q", 0 );
  int readInsertFile = P.getOptionIntValue( "-i", 1 );

  int LEAVE_WRAP = 32;
  parlay::sequence<PointType<coord, 15>> wp;
  // parlay::sequence<PointID<coord, 15>> wp;
  std::string name, insertFile = "";

  //* initialize points
  if ( iFile != NULL ) {
    name = std::string( iFile );
    name = name.substr( name.rfind( "/" ) + 1 );
    std::cout << name << " ";
    auto [n, d] = read_points<PointType<coord, 15>>( iFile, wp, K );
    // auto [n, d] = read_points<PointID<coord, 15>>( iFile, wp, K );
    N = n;
    assert( d == Dim );
  } else {  //* construct data byself
    K = 100;
    generate_random_points<PointType<coord, 15>>( wp, 1000000, N, Dim );
    // generate_random_points<PointID<coord, 15>>( wp, 1000000, N, Dim );
    std::string name = std::to_string( N ) + "_" + std::to_string( Dim ) + ".in";
    std::cout << name << " ";
  }

  if ( readInsertFile == 1 ) {
    int id = std::stoi( name.substr( 0, name.find_first_of( '.' ) ) );
    id = ( id + 1 ) % 3;  //! MOD graph number used to test
    if ( !id ) id++;
    int pos = std::string( iFile ).rfind( "/" ) + 1;
    insertFile = std::string( iFile ).substr( 0, pos ) + std::to_string( id ) + ".in";
  }

  assert( N > 0 && Dim > 0 && K > 0 && LEAVE_WRAP >= 1 );

  if ( tag == -1 ) {
    //* serial run
    // todo rewrite test serial code
    // testSerialKDtree( Dim, LEAVE_WRAP, wp, N, K );
  } else if ( Dim == 2 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 2> {
      return PointType<coord, 2>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 2>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile, tag,
                                      queryType );
  } else if ( Dim == 3 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 3> {
      return PointType<coord, 3>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 3>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile, tag,
                                      queryType );
  } else if ( Dim == 5 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 5> {
      return PointType<coord, 5>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 5>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile, tag,
                                      queryType );
  } else if ( Dim == 7 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 7> {
      return PointType<coord, 7>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 7>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile, tag,
                                      queryType );
  } else if ( Dim == 9 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 9> {
      return PointType<coord, 9>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 9>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile, tag,
                                      queryType );
  } else if ( Dim == 10 ) {
    auto pts = parlay::tabulate( N, [&]( size_t i ) -> PointType<coord, 10> {
      return PointType<coord, 10>( wp[i].pnt.begin() );
    } );
    decltype( wp )().swap( wp );
    testOctTree<PointType<coord, 10>>( Dim, LEAVE_WRAP, pts, N, K, rounds, insertFile,
                                       tag, queryType );
  }

  return 0;
}
