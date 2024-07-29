#pragma once

#include "../base_tree.h"

namespace cpdd {

template<typename Point, uint8_t kBDO>
inline size_t BaseTree<Point, kBDO>::GetImbalanceRatio() {
    if (const auto env_p = std::getenv("kInbalanceRatio")) {
        return static_cast<size_t>(std::stoi(env_p));
    } else {
        return static_cast<size_t>(kInbalanceRatio);
    }
}

template<typename Point, uint8_t kBDO>
inline bool BaseTree<Point, kBDO>::ImbalanceNode(const size_t l,
                                                 const size_t n) {
    if (n == 0) return true;
    return Num::Gt(static_cast<size_t>(std::abs(100.0 * l / n - 50.0)),
                   GetImbalanceRatio());
}

// template<typename Point, uint8_t kBDO>
// template<typename SplitRule>
// inline BaseTree<Point, kBDO>::DimsType BaseTree<Point,
// kBDO>::PickRebuildDim(
//     const Node* T, const DimsType d, const DimsType DIM) {
//     PickRebuildDim_(T, d, DIM, SplitRule());
// }
//
// template<typename Point, uint8_t kBDO>
// inline BaseTree<Point, kBDO>::DimsType BaseTree<Point,
// kBDO>::PickRebuildDim_(
//     const Node* T, const DimsType d, const DimsType DIM, MaxStretchDimTag) {
//     return 0;
// }
//
// template<typename Point, uint8_t kBDO>
// inline BaseTree<Point, kBDO>::DimsType BaseTree<Point,
// kBDO>::PickRebuildDim_(
//     const Node* T, const DimsType d, const DimsType DIM, RotateDimTag) {
//     return d;
// }

// template<typename Point, uint8_t kBDO>
// inline BaseTree<Point, kBDO>::DimsType BaseTree<Point,
// kBDO>::PickMaxStretchDim(
//     const Box& bx) {
//     DimsType d(0);
//     Coord diff(bx.second.pnt[0] - bx.first.pnt[0]);
//     assert(Num::Geq(diff, 0));
//     for (int i = 1; i < bx.first.get_dim(); i++) {
//         if (Num::Gt(bx.second.pnt[i] - bx.first.pnt[i], diff)) {
//             diff = bx.second.pnt[i] - bx.first.pnt[i];
//             d = i;
//         }
//     }
//     return d;
// }
}  // namespace cpdd
