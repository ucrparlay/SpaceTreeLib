#include "../base_tree.h"

namespace cpdd {
// DimsType d =
//     (split_rule_ == kMaxStretchDim ? pick_max_stretch_dim(bx, DIM) :
//     dim);
//
template<typename Point>
struct BaseSplitterRule {
    using BT = BaseTree<Point>;
    using Slice = BT::Slice;
    using DimsType = BT::DimsType;
    using Box = BT::Box;
    using Num = BT::Num;
    using Coord = BT::Coord;

    virtual DimsType FindCuttingDimension(const Box& bx, const DimsType dim,
                                          const DimsType DIM) = 0;

    virtual std::pair<Box, DimsType> SwitchDimension(const Slice In,
                                                     const DimsType dim,
                                                     const DimsType DIM,
                                                     const Box& bx) = 0;

    virtual DimsType FindRebuildDimension(const DimsType dim,
                                          const DimsType DIM) = 0;
};

template<typename Point>
struct MaxStretchDim : BaseSplitterRule<Point> {
    using BSR = BaseSplitterRule<Point>;
    using BT = BSR::BT;
    using Slice = BSR::Slice;
    using DimsType = BSR::DimsType;
    using Box = BSR::Box;
    using Coord = BSR::Coord;
    using Num = BSR::Num;

    DimsType FindCuttingDimension(const Box& bx,
                                  [[maybe_unused]] const DimsType dim,
                                  const DimsType DIM) override {
        DimsType d(0);
        Coord diff(bx.second.pnt[0] - bx.first.pnt[0]);
        assert(Num::Geq(diff, 0));
        for (DimsType i = 1; i < DIM; ++i) {
            if (Num::Gt(bx.second.pnt[i] - bx.first.pnt[i], diff)) {
                diff = bx.second.pnt[i] - bx.first.pnt[i];
                d = i;
            }
        }
        return d;
    };

    DimsType FindRebuildDimension(
        [[maybe_unused]] const DimsType dim,
        [[maybe_unused]] const DimsType DIM) override {
        return 0;
    };

    std::pair<Box, DimsType> SwitchDimension(
        const Slice In, const DimsType dim, [[maybe_unused]] const DimsType DIM,
        [[maybe_unused]] const Box& bx) override {
        return std::make_pair(BT::GetBox(In), dim);
    };
};

template<typename Point>
struct RotateDim : BaseSplitterRule<Point> {
    using BSR = BaseSplitterRule<Point>;
    using BT = BSR::BT;
    using Slice = BSR::Slice;
    using DimsType = BSR::DimsType;
    using Box = BSR::Box;
    using Coord = BSR::Coord;
    using Num = BSR::Num;
    using PointsIter = BT::PointsIter;

    DimsType FindCuttingDimension(
        [[maybe_unused]] const Box& bx, const DimsType dim,
        [[maybe_unused]] const DimsType DIM) override {
        return dim;
    };

    DimsType FindRebuildDimension(
        [[maybe_unused]] const DimsType dim,
        [[maybe_unused]] const DimsType DIM) override {
        return dim;
    };

    std::pair<Box, DimsType> SwitchDimension(
        const Slice In, [[maybe_unused]] const DimsType dim, const DimsType DIM,
        [[maybe_unused]] const Box& bx) override {
        DimsType d(0);
        for (DimsType i = 0; i < DIM; ++i) {
            if (!Num::Eq(In.begin()->pnt[i], (--In.end())->pnt[i])) {
                d = i;
            }
        }
        return std::make_pair(bx, d);
    };
};
}  // namespace cpdd
