#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point, typename DerivedTree, uint_fast8_t kBDO>
struct BaseTree<Point, DerivedTree, kBDO>::BoxCut {
    using BT = BaseTree<Point, DerivedTree, kBDO>;

    BoxCut(const Box& box, const HyperPlane& hp, bool go_left) :
        box(box), hp(hp), go_left(go_left) {}

    inline const Box& GetFirstBoxCut() {
        mod_dim =
            go_left ? &box.second.pnt[hp.second] : &box.first.pnt[hp.second];
        std::ranges::swap(hp.first, *mod_dim);
        return box;
    }

    inline const Box& GetSecondBoxCut() {
        std::ranges::swap(hp.first, *mod_dim);
        mod_dim =
            go_left ? &box.first.pnt[hp.second] : &box.second.pnt[hp.second];
        *mod_dim = hp.first;
        return box;
    }

    inline const Box& GetBox() const { return box; }

    Box box;
    Coord* mod_dim;
    HyperPlane hp;  // PARA: the split and the cutting dimension
    const bool go_left;
};
};  // namespace cpdd
