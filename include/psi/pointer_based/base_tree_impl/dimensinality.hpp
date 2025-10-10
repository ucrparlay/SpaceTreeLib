#ifndef PSI_POINTER_BASED_BASE_TREE_IMPL_DIMENSINALITY_HPP_
#define PSI_POINTER_BASED_BASE_TREE_IMPL_DIMENSINALITY_HPP_

#include "../base_tree.h"

namespace psi {
namespace pointer_based {

template <class TypeTrait, typename DerivedTree>
inline size_t BaseTree<TypeTrait, DerivedTree>::GetImbalanceRatio() {
  return static_cast<size_t>(kInbalanceRatio);
  // if (auto const env_p = std::getenv("kInbalanceRatio")) {
  //   return static_cast<size_t>(std::stoi(env_p));
  // } else {
  //   return static_cast<size_t>(kInbalanceRatio);
  // }
}

template <class TypeTrait, typename DerivedTree>
inline bool BaseTree<TypeTrait, DerivedTree>::ImbalanceNode(size_t const l,
                                                            size_t const n) {
  if (n == 0) return true;
  return Num::Gt(
      static_cast<size_t>(std::abs(
          100.0 * static_cast<double>(l) / static_cast<double>(n) - 50.0)),
      GetImbalanceRatio());
}

template <class TypeTrait, typename DerivedTree>
inline bool BaseTree<TypeTrait, DerivedTree>::SparcyNode(size_t const rm,
                                                         size_t const n) {
  // PERF: to avoid the case that the new leaf is about 32 and then next slight
  // larger insert will break the leaf
  return n - rm < kThinLeaveWrap;
}
}  // namespace pointer_based
}  // namespace psi

#endif  // PSI_POINTER_BASED_BASE_TREE_IMPL_DIMENSINALITY_HPP_