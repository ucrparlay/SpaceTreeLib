#pragma once

#include "../base_tree.h"

namespace cpdd {

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline size_t BaseTree<Point, DerivedTree, kBDO>::GetImbalanceRatio() {
    if (const auto env_p = std::getenv("kInbalanceRatio")) {
        return static_cast<size_t>(std::stoi(env_p));
    } else {
        return static_cast<size_t>(kInbalanceRatio);
    }
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::ImbalanceNode(const size_t l,
                                                              const size_t n) {
    if (n == 0) return true;
    return Num::Gt(static_cast<size_t>(std::abs(100.0 * l / n - 50.0)),
                   GetImbalanceRatio());
}

template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::SparcyNode(const size_t rm,
                                                           const size_t n) {
    return n - rm < kLeaveWrap;
}
}  // namespace cpdd
