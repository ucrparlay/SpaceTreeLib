#pragma once

#include "../baseTree.h"

namespace cpdd {

template<typename point>
template<typename Slice>
void
baseTree<point>::flatten( typename baseTree<point>::node* T, Slice Out,
                          bool granularity ) {
  assert( T->size == Out.size() );
  if ( T->size == 0 ) return;

  if ( T->is_leaf ) {
    leaf* TL = static_cast<leaf*>( T );
    for ( int i = 0; i < TL->size; i++ ) {
      Out[i] = TL->pts[( !T->is_dummy ) * i];
    }
    return;
  }

  interior* TI = static_cast<interior*>( T );
  assert( TI->size == TI->left->size + TI->right->size );
  parlay::par_do_if(
      ( granularity && TI->size > SERIAL_BUILD_CUTOFF ) || !granularity,
      [&]() { flatten( TI->left, Out.cut( 0, TI->left->size ), granularity ); },
      [&]() { flatten( TI->right, Out.cut( TI->left->size, TI->size ), granularity ); } );

  return;
}

template<typename point>
inline void
baseTree<point>::update_interior( typename baseTree<point>::node* T,
                                  typename baseTree<point>::node* L,
                                  typename baseTree<point>::node* R ) {
  assert( !T->is_leaf );
  interior* TI = static_cast<interior*>( T );
  TI->size = L->size + R->size;
  TI->left = L;
  TI->right = R;
  return;
}

template<typename point>
uint_fast8_t
baseTree<point>::retrive_tag( const point& p, const node_tags& tags ) {
  uint_fast8_t k = 1;
  interior* TI;
  while ( k <= PIVOT_NUM && ( !tags[k].first->is_leaf ) ) {
    TI = static_cast<interior*>( tags[k].first );
    k = Num::Lt( p.pnt[TI->split.second], TI->split.first ) ? k << 1 : k << 1 | 1;
  }
  assert( tags[k].second < BUCKET_NUM );
  return tags[k].second;
}

template<typename point>
void
baseTree<point>::seieve_points( slice A, slice B, const size_t n, const node_tags& tags,
                                parlay::sequence<balls_type>& sums,
                                const bucket_type tagsNum ) {
  size_t num_block = ( n + BLOCK_SIZE - 1 ) >> LOG2_BASE;
  parlay::sequence<parlay::sequence<balls_type>> offset(
      num_block, parlay::sequence<balls_type>( tagsNum ) );
  assert( offset.size() == num_block && offset[0].size() == tagsNum &&
          offset[0][0] == 0 );
  parlay::parallel_for( 0, num_block, [&]( size_t i ) {
    for ( size_t j = i << LOG2_BASE; j < std::min( ( i + 1 ) << LOG2_BASE, n ); j++ ) {
      offset[i][std::move( retrive_tag( A[j], tags ) )]++;
    }
  } );

  sums = parlay::sequence<balls_type>( tagsNum );
  for ( size_t i = 0; i < num_block; i++ ) {
    auto t = offset[i];
    offset[i] = sums;
    for ( int j = 0; j < tagsNum; j++ ) {
      sums[j] += t[j];
    }
  }

  parlay::parallel_for( 0, num_block, [&]( size_t i ) {
    auto v = parlay::sequence<balls_type>::uninitialized( tagsNum );
    int tot = 0, s_offset = 0;
    for ( int k = 0; k < tagsNum - 1; k++ ) {
      v[k] = tot + offset[i][k];
      tot += sums[k];
      s_offset += offset[i][k];
    }
    v[tagsNum - 1] = tot + ( ( i << LOG2_BASE ) - s_offset );
    for ( size_t j = i << LOG2_BASE; j < std::min( ( i + 1 ) << LOG2_BASE, n ); j++ ) {
      B[v[std::move( retrive_tag( A[j], tags ) )]++] = A[j];
    }
  } );

  return;
}

template<typename point>
typename baseTree<point>::node*
baseTree<point>::delete_tree() {
  if ( this->root == nullptr ) {
    return this->root;
  }
  delete_tree_recursive( this->root );
  this->root = nullptr;
  return this->root;
}

template<typename point>  //* delete tree in parallel
void
baseTree<point>::delete_tree_recursive( node* T, bool granularity ) {
  if ( T == nullptr ) return;
  if ( T->is_leaf ) {
    free_node<point, leaf>( T );
  } else {
    interior* TI = static_cast<interior*>( T );

    // NOTE: enable granularity control by default, if it is disabled, always delete in
    // parallel
    parlay::par_do_if( ( granularity && T->size > SERIAL_BUILD_CUTOFF ) || !granularity,
                       [&] { delete_tree_recursive( TI->left, granularity ); },
                       [&] { delete_tree_recursive( TI->right, granularity ); } );
    free_node<point, interior>( T );
  }
}

}  // namespace cpdd
