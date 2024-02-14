#pragma once

#include "../baseTree.h"

namespace cpdd {
template<typename point>
inline baseTree<point>::dim_type
baseTree<point>::pick_rebuild_dim( const node* T, const dim_type d, const dim_type DIM ) {
  if ( this->_split_rule == MAX_STRETCH_DIM ) {
    return 0;
  } else if ( this->_split_rule == ROTATE_DIM ) {
    return d;
  } else {
    // WARN: customize it before using
    return 0;
  }
}

template<typename point>
inline baseTree<point>::dim_type
baseTree<point>::pick_max_stretch_dim( const box& bx, const dim_type DIM ) {
  dim_type d( 0 );
  coord diff( bx.second.pnt[0] - bx.first.pnt[0] );
  assert( Num::Geq( diff, 0 ) );
  for ( int i = 1; i < DIM; i++ ) {
    if ( Num::Gt( bx.second.pnt[i] - bx.first.pnt[i], diff ) ) {
      diff = bx.second.pnt[i] - bx.first.pnt[i];
      d = i;
    }
  }
  return d;
}
}  // namespace cpdd
