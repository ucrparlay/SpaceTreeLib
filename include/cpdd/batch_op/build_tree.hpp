#pragma once

#include <cassert>
#include "../baseTree.h"
#include "parlay/primitives.h"

namespace cpdd {

// template<typename point>
// void
// baseTree<point>::build( slice A, const dim_type DIM ) {
//   points B = points::uninitialized( A.size() );
//   this->bbox = get_box( A );
//   this->root = build_recursive( A, B.cut( 0, A.size() ), 0, DIM, this->bbox );
//   assert( this->root != nullptr );
//   return;
// }
//

template<typename point>
void
baseTree<point>::divide_rotate( slice In, splitter_s& pivots, dim_type dim,
                                bucket_type idx, bucket_type deep, bucket_type& bucket,
                                const dim_type DIM, box_s& boxs, const box& bx ) {
  if ( deep > BUILD_DEPTH_ONCE ) {
    //! sometimes cut dimension can be -1
    //! never use pivots[idx].first to check whether it is in bucket; instead,
    //! use idx > PIVOT_NUM
    boxs[bucket] = bx;
    pivots[idx] = splitter( -1, bucket++ );
    return;
  }
  size_t n = In.size();
  uint_fast8_t d =
      ( _split_rule == MAX_STRETCH_DIM ? pick_max_stretch_dim( bx, DIM ) : dim );
  assert( d >= 0 && d < DIM );

  std::ranges::nth_element( In.begin(), In.begin() + n / 2, In.end(),
                            [&]( const point& p1, const point& p2 ) {
                              return Num::Lt( p1.pnt[d], p2.pnt[d] );
                            } );
  pivots[idx] = splitter( In[n / 2].pnt[d], d );

  box lbox( bx ), rbox( bx );
  lbox.second.pnt[d] = pivots[idx].first;  //* loose
  rbox.first.pnt[d] = pivots[idx].first;

  d = ( d + 1 ) % DIM;
  divide_rotate( In.cut( 0, n / 2 ), pivots, d, 2 * idx, deep + 1, bucket, DIM, boxs,
                 lbox );
  divide_rotate( In.cut( n / 2, n ), pivots, d, 2 * idx + 1, deep + 1, bucket, DIM, boxs,
                 rbox );
  return;
}

//* starting at dimesion dim and pick pivots in a rotation manner
template<typename point>
void
baseTree<point>::pick_pivots( slice In, const size_t& n, splitter_s& pivots,
                              const dim_type dim, const dim_type DIM, box_s& boxs,
                              const box& bx ) {
  size_t size = std::min( n, (size_t)32 * BUCKET_NUM );
  assert( size <= n );

  //* samples
  points arr = points::uninitialized( size );
  auto indexs = parlay::sequence<uint64_t>::uninitialized( size );
  for ( size_t i = 0; i < size; i++ ) {
    indexs[i] = _hash64( i ) % n;
  }
  std::sort( indexs.begin(), indexs.end() );
  for ( size_t i = 0; i < size; i++ ) {
    arr[i] = In[indexs[i]];
  }

  //* pick pivots
  bucket_type bucket = 0;
  divide_rotate( arr.cut( 0, size ), pivots, dim, 1, 1, bucket, DIM, boxs, bx );
  assert( bucket == BUCKET_NUM );
  return;
}

template<typename point>
inline uint_fast8_t
baseTree<point>::find_bucket( const point& p, const splitter_s& pivots ) {
  uint_fast8_t k = 1;
  while ( k <= PIVOT_NUM ) {
    // k = p.pnt[pivots[k].second] < pivots[k].first ? k << 1 : k << 1 | 1;
    k = Num::Lt( p.pnt[pivots[k].second], pivots[k].first ) ? k << 1 : k << 1 | 1;
  }
  assert( pivots[k].first == -1 );
  return pivots[k].second;
}

template<typename point>
void
baseTree<point>::partition( slice A, slice B, const size_t n, const splitter_s& pivots,
                            parlay::sequence<balls_type>& sums ) {
  size_t num_block = ( n + BLOCK_SIZE - 1 ) >> LOG2_BASE;
  parlay::sequence<parlay::sequence<balls_type>> offset(
      num_block, parlay::sequence<balls_type>( BUCKET_NUM ) );
  assert( offset.size() == num_block && offset[0].size() == BUCKET_NUM &&
          offset[0][0] == 0 );
  parlay::parallel_for( 0, num_block, [&]( size_t i ) {
    for ( size_t j = i << LOG2_BASE; j < std::min( ( i + 1 ) << LOG2_BASE, n ); j++ ) {
      offset[i][std::move( find_bucket( A[j], pivots ) )]++;
    }
  } );

  sums = parlay::sequence<balls_type>( BUCKET_NUM );
  for ( size_t i = 0; i < num_block; i++ ) {
    auto t = offset[i];
    offset[i] = sums;
    for ( bucket_type j = 0; j < BUCKET_NUM; j++ ) {
      sums[j] += t[j];
    }
  }

  parlay::parallel_for( 0, num_block, [&]( size_t i ) {
    auto v = parlay::sequence<balls_type>::uninitialized( BUCKET_NUM );
    size_t tot = 0, s_offset = 0;
    for ( bucket_type k = 0; k < BUCKET_NUM - 1; k++ ) {
      v[k] = tot + offset[i][k];
      tot += sums[k];
      s_offset += offset[i][k];
    }
    v[BUCKET_NUM - 1] = tot + ( ( i << LOG2_BASE ) - s_offset );
    for ( size_t j = i << LOG2_BASE; j < std::min( ( i + 1 ) << LOG2_BASE, n ); j++ ) {
      uint_fast8_t k = find_bucket( A[j], pivots );
      // B[v[std::move( find_bucket( A[j], pivots ) )]++] = A[j];
      B[v[k]++] = A[j];
    }
  } );

  return;
}

// template<typename point>
// typename baseTree<point>::node*
// baseTree<point>::build_inner_tree( bucket_type idx, splitter_s& pivots,
//                                    parlay::sequence<node*>& treeNodes ) {
//   if ( idx > PIVOT_NUM ) {
//     assert( idx - PIVOT_NUM - 1 < BUCKET_NUM );
//     return treeNodes[idx - PIVOT_NUM - 1];
//   }
//   node *L, *R;
//   L = build_inner_tree( idx << 1, pivots, treeNodes );
//   R = build_inner_tree( idx << 1 | 1, pivots, treeNodes );
//   return alloc_interior_node( L, R, pivots[idx] );
// }

template<typename point>
typename baseTree<point>::points_iter
baseTree<point>::serial_partition( slice In, dim_type d ) {
  size_t n = In.size();
  std::ranges::nth_element( In.begin(), In.begin() + n / 2, In.end(),
                            [&]( const point& p1, const point& p2 ) {
                              return Num::Lt( p1.pnt[d], p2.pnt[d] );
                            } );
  std::ranges::subrange _2ndGroup = std::ranges::partition(
      In.begin(), In.begin() + n / 2,
      [&]( const point& p ) { return Num::Lt( p.pnt[d], In[n / 2].pnt[d] ); } );
  if ( _2ndGroup.begin() == In.begin() ) {  //* handle duplicated medians
    _2ndGroup = std::ranges::partition(
        In.begin() + n / 2, In.end(),
        [&]( const point& p ) { return Num::Eq( p.pnt[d], In[n / 2].pnt[d] ); } );
  }
  return _2ndGroup.begin();
}

// template<typename point>
// typename baseTree<point>::node*
// baseTree<point>::serial_build_recursive( slice In, slice Out, dim_type dim,
//                                          const dim_type DIM, const box& bx ) {
//   size_t n = In.size();
//   if ( n == 0 ) return alloc_empty_leaf();
//   if ( n <= LEAVE_WRAP ) return alloc_leaf_node( In );
//
//   auto parallel_partition = [&]( slice In, dim_type d ) {
//     auto cmp = [&]( const point& p1, const point& p2 ) {
//       return Num::Lt( p1.pnt[d], p2.pnt[d] );
//     };
//     auto iter_pivot = parlay::kth_smallest( In, n / 2, cmp );
//     auto iter_min = parlay::kth_smallest( In, 0, cmp );
//     auto ucmp = *iter_pivot == *iter_min ? &Num::Leq : &Num::Lt;
//     auto res = parlay::internal::split_two(
//         In, parlay::delayed_seq<bool>( n, [&]( const point& p ) {
//           return ( *ucmp )( p.pnt[d], iter_pivot->pnt[d] );
//         } ) );
//     parlay::copy( std::move( res.first ), In );
//     return In.begin() + res.second;
//   };
//
//   dim_type d = ( _split_rule == MAX_STRETCH_DIM ? pick_max_stretch_dim( bx, DIM ) : dim
//   );
//   // points_iter splitIter =
//   //     n > 5000000 ? parallel_partition( In, d ) : serial_partition( In, d );
//   points_iter splitIter = serial_partition( In, d );
//   points_iter diffEleIter;
//
//   splitter split;
//
//   if ( splitIter <= In.begin() + n / 2 ) {
//     split = splitter( In[n / 2].pnt[d], d );
//   } else if ( splitIter != In.end() ) {
//     auto minEleIter = std::ranges::min_element( splitIter, In.end(),
//                                                 [&]( const point& p1, const point& p2 )
//                                                 {
//                                                   return Num::Lt( p1.pnt[d], p2.pnt[d]
//                                                   );
//                                                 } );
//     split = splitter( minEleIter->pnt[d], d );
//   } else if ( In.end() ==
//               ( diffEleIter = std::ranges::find_if_not( In, [&]( const point& p ) {
//                   return p.sameDimension( In[0] );
//                 } ) ) ) {  //* check whether all elements are identical
//     return alloc_dummy_leaf( In );
//   } else {  //* current dim d is same but other dims are not
//     if ( _split_rule == MAX_STRETCH_DIM ) {  //* next recursion redirects to new dim
//       return serial_build_recursive( In, Out, d, DIM, get_box( In ) );
//     } else if ( _split_rule == ROTATE_DIM ) {  //* switch dim, break rotation order
//       points_iter compIter =
//           diffEleIter == In.begin() ? std::ranges::prev( In.end() ) : In.begin();
//       assert( compIter != diffEleIter );
//       for ( int i = 0; i < DIM; i++ ) {
//         if ( !Num::Eq( diffEleIter->pnt[i], compIter->pnt[i] ) ) {
//           d = i;
//           break;
//         }
//       }
//       return serial_build_recursive( In, Out, d, DIM, bx );
//     }
//   }
//
//   assert( std::ranges::all_of(
//               In.begin(), splitIter,
//               [&]( point& p ) { return Num::Lt( p.pnt[split.second], split.first ); } )
//               &&
//           std::ranges::all_of( splitIter, In.end(), [&]( point& p ) {
//             return Num::Geq( p.pnt[split.second], split.first );
//           } ) );
//
//   box lbox( bx ), rbox( bx );
//   lbox.second.pnt[d] = split.first;  //* loose
//   rbox.first.pnt[d] = split.first;
//
//   d = ( d + 1 ) % DIM;
//   node *L, *R;
//
//   L = serial_build_recursive( In.cut( 0, splitIter - In.begin() ),
//                               Out.cut( 0, splitIter - In.begin() ), d, DIM, lbox );
//   R = serial_build_recursive( In.cut( splitIter - In.begin(), n ),
//                               Out.cut( splitIter - In.begin(), n ), d, DIM, rbox );
//   return alloc_interior_node( L, R, split );
// }

// template<typename point>
// typename baseTree<point>::node*
// baseTree<point>::build_recursive( slice In, slice Out, dim_type dim, const dim_type
// DIM,
//                                   const box& bx ) {
//   assert( In.size() == 0 || within_box( get_box( In ), bx ) );
//
//   // if ( In.size() ) {
//   if ( In.size() <= SERIAL_BUILD_CUTOFF ) {
//     return serial_build_recursive( In, Out, dim, DIM, bx );
//   }
//
//   //* parallel partitons
//   auto pivots = splitter_s::uninitialized( PIVOT_NUM + BUCKET_NUM + 1 );
//   auto boxs = box_s::uninitialized( BUCKET_NUM );
//   parlay::sequence<balls_type> sums;
//
//   pick_pivots( In, In.size(), pivots, dim, DIM, boxs, bx );
//   partition( In, Out, In.size(), pivots, sums );
//
//   //* if random sampling failed to split points, re-partitions using serail approach
//   auto treeNodes = parlay::sequence<node*>::uninitialized( BUCKET_NUM );
//   auto nodesMap = parlay::sequence<bucket_type>::uninitialized( BUCKET_NUM );
//
//   bucket_type zeros = 0, cnt = 0;
//   for ( bucket_type i = 0; i < BUCKET_NUM; i++ ) {
//     if ( !sums[i] ) {
//       zeros++;
//       treeNodes[i] = alloc_empty_leaf();
//     } else {
//       nodesMap[cnt++] = i;
//     }
//   }
//
//   if ( zeros == BUCKET_NUM - 1 ) {  // * switch to seral
//                                     // TODO: add parallelsim within this call
//     return serial_build_recursive( In, Out, dim, DIM, bx );
//   }
//
//   dim = ( dim + BUILD_DEPTH_ONCE ) % DIM;
//
//   parlay::parallel_for(
//       0, BUCKET_NUM - zeros,
//       [&]( bucket_type i ) {
//         size_t start = 0;
//         for ( bucket_type j = 0; j < nodesMap[i]; j++ ) {
//           start += sums[j];
//         }
//
//         treeNodes[nodesMap[i]] = build_recursive(
//             Out.cut( start, start + sums[nodesMap[i]] ),
//             In.cut( start, start + sums[nodesMap[i]] ), dim, DIM, boxs[nodesMap[i]] );
//       },
//       1 );
//   return build_inner_tree( 1, pivots, treeNodes );
// }

}  // namespace cpdd
