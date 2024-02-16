#pragma once

#include "../base_tree.h"

namespace cpdd {

template<typename point>
inline size_t
baseTree<point>::get_imbalance_ratio() {
  if (const auto env_p = std::getenv("INBALANCE_RATIO")) {
    return static_cast<size_t>(std::stoi(env_p));
  } else {
    return static_cast<size_t>(INBALANCE_RATIO);
  }
}

template<typename point>
inline bool
baseTree<point>::inbalance_node(const size_t l, const size_t n) {
  if (n == 0) return true;
  return Num::Gt(static_cast<size_t>(std::abs(100.0 * l / n - 50.0)),
                 get_imbalance_ratio());
}

template<typename point>
inline baseTree<point>::dim_type
baseTree<point>::pick_rebuild_dim(const node* T, const dim_type d,
                                  const dim_type DIM) {
  if (this->_split_rule == MAX_STRETCH_DIM) {
    return 0;
  } else if (this->_split_rule == ROTATE_DIM) {
    return d;
  } else {
    // WARN: customize it before using
    return 0;
  }
}

template<typename point>
inline baseTree<point>::dim_type
baseTree<point>::pick_max_stretch_dim(const box& bx, const dim_type DIM) {
  dim_type d(0);
  coord diff(bx.second.pnt[0] - bx.first.pnt[0]);
  assert(Num::Geq(diff, 0));
  for (int i = 1; i < DIM; i++) {
    if (Num::Gt(bx.second.pnt[i] - bx.first.pnt[i], diff)) {
      diff = bx.second.pnt[i] - bx.first.pnt[i];
      d = i;
    }
  }
  return d;
}
}  // namespace cpdd
