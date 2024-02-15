#pragma once

#include <algorithm>
#include "../baseTree.h"

namespace cpdd {

template<typename point>
template<typename StoreType>
size_t
baseTree<point>::range_query_parallel( const typename baseTree<point>::box& queryBox,
                                       StoreType Out, double& tim ) {
  parlay::internal::timer t;
  tim = 0.0;
  t.start();
  auto ST = range_count_save_path( this->root, queryBox, this->bbox );
  size_t sz = ST->size;
  tim += t.next_time();

  range_query_recursive_parallel( this->root, ST, Out.cut( 0, sz ), queryBox );

  t.start();
  delete_simple_tree_recursive( ST );
  t.stop();
  tim += t.next_time();
  return sz;
}

template<typename point>
template<typename StoreType>
size_t
baseTree<point>::range_query_serial( const typename baseTree<point>::box& queryBox,
                                     StoreType Out ) {
  size_t s = 0;
  range_query_recursive_serial( this->root, Out, s, queryBox, this->bbox );
  return s;
}

template<typename point>
template<typename StoreType>
void
baseTree<point>::range_query_recursive_serial( node* T, StoreType Out, size_t& s,
                                               const box& queryBox, const box& nodeBox ) {
  if ( T->is_leaf ) {
    leaf* TL = static_cast<leaf*>( T );
    if ( T->is_dummy ) {
      if ( within_box( TL->pts[0], queryBox ) ) {
        for ( int i = 0; i < TL->size; i++ )
          Out[s++] = TL->pts[0];
      }
    } else {
      for ( int i = 0; i < TL->size; i++ )
        if ( within_box( TL->pts[i], queryBox ) ) {
          Out[s++] = TL->pts[i];
        }
    }
    return;
  }

  interior* TI = static_cast<interior*>( T );
  // box lbox( nodeBox ), rbox( nodeBox );
  box abox( nodeBox );
  // lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
  // rbox.first.pnt[TI->split.second] = TI->split.first;

  auto recurse = [&]( node* Ts, const box& bx ) -> void {
    if ( !box_intersect_box( bx, queryBox ) ) {
      return;
    } else if ( within_box( bx, queryBox ) ) {
      flatten( Ts, Out.cut( s, s + Ts->size ) );
      s += Ts->size;
      return;
    } else {
      range_query_recursive_serial( Ts, Out, s, queryBox, bx );
      return;
    }
  };

  auto& mod_dim = abox.second.pnt[TI->split.second];
  auto split = TI->split.first;
  std::swap( mod_dim, split );
  recurse( TI->left, abox );

  std::swap( mod_dim, split );
  abox.first.pnt[TI->split.second] = split;
  recurse( TI->right, abox );

  return;
}

}  // namespace cpdd
