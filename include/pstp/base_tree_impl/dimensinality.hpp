#ifndef PSTP_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
#define PSTP_BASE_TREE_IMPL_DIMENSIONALITY_HPP_

#include "../base_tree.h"

namespace pstp {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline size_t
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetImbalanceRatio() {
  return static_cast<size_t>(kInbalanceRatio);
  // if (auto const env_p = std::getenv("kInbalanceRatio")) {
  //   return static_cast<size_t>(std::stoi(env_p));
  // } else {
  //   return static_cast<size_t>(kInbalanceRatio);
  // }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::ImbalanceNode(
    size_t const l, size_t const n) {
  if (n == 0) return true;
  return Num::Gt(
      static_cast<size_t>(std::abs(
          100.0 * static_cast<double>(l) / static_cast<double>(n) - 50.0)),
      GetImbalanceRatio());
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::SparcyNode(
    size_t const rm, size_t const n) {
  // PERF: to avoid the case that the new leaf is about 32 and then next slight
  // larger insert will break the leaf
  return n - rm < kThinLeaveWrap;
}
}  // namespace pstp

#endif  // PSTP_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
