#pragma once

#include "../base_tree.h"

namespace cpdd {

template<typename Point>
inline size_t BaseTree<Point>::GetImbalanceRatio() {
    if (const auto env_p = std::getenv("kInbalanceRatio")) {
        return static_cast<size_t>(std::stoi(env_p));
    } else {
        return static_cast<size_t>(kInbalanceRatio);
    }
}

template<typename Point>
inline bool BaseTree<Point>::ImbalanceNode(const size_t l, const size_t n) {
    if (n == 0) return true;
    return Num::Gt(static_cast<size_t>(std::abs(100.0 * l / n - 50.0)), GetImbalanceRatio());
}

template<typename Point>
template<typename SplitRule>
inline BaseTree<Point>::DimsType BaseTree<Point>::PickRebuildDim(const Node* T, const DimsType d, const DimsType DIM) {
    PickRebuildDim_(T, d, DIM, SplitRule());
}

template<typename Point>
inline BaseTree<Point>::DimsType BaseTree<Point>::PickRebuildDim_(const Node* T, const DimsType d, const DimsType DIM,
                                                                  MaxStretchDimTag) {
    return 0;
}

template<typename Point>
inline BaseTree<Point>::DimsType BaseTree<Point>::PickRebuildDim_(const Node* T, const DimsType d, const DimsType DIM,
                                                                  RotateDimTag) {
    return d;
}

template<typename Point>
inline BaseTree<Point>::DimsType BaseTree<Point>::PickMaxStretchDim(const Box& bx, const DimsType DIM) {
    DimsType d(0);
    Coord diff(bx.second.pnt[0] - bx.first.pnt[0]);
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
