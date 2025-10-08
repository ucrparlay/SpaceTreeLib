#ifndef PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
#define PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_

#include "../base_tree.h"

#define BASETREE_TEMPLATE                                                 \
  template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define BASETREE_CLASS BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>

namespace psi {

BASETREE_TEMPLATE
inline size_t BASETREE_CLASS::GetImbalanceRatio() {
  return static_cast<size_t>(kInbalanceRatio);
  // if (auto const env_p = std::getenv("kInbalanceRatio")) {
  //   return static_cast<size_t>(std::stoi(env_p));
  // } else {
  //   return static_cast<size_t>(kInbalanceRatio);
  // }
}

BASETREE_TEMPLATE
inline bool BASETREE_CLASS::ImbalanceNode(size_t const l, size_t const n) {
  if (n == 0) return true;
  return Num::Gt(
      static_cast<size_t>(std::abs(
          100.0 * static_cast<double>(l) / static_cast<double>(n) - 50.0)),
      GetImbalanceRatio());
}

BASETREE_TEMPLATE
inline bool BASETREE_CLASS::SparcyNode(size_t const rm, size_t const n) {
  // PERF: to avoid the case that the new leaf is about 32 and then next slight
  // larger insert will break the leaf
  return n - rm < kThinLeaveWrap;
}
}  // namespace psi

#undef BASETREE_TEMPLATE
#undef BASETREE_CLASS

#endif  // PSI_BASE_TREE_IMPL_DIMENSIONALITY_HPP_
