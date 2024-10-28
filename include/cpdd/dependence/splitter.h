#pragma once

#include "../base_tree.h"

namespace cpdd {
template <typename Point>
struct BaseSplitterRule {
  using BT = BaseTree<Point, BaseSplitterRule<Point>>;
  using Slice = BT::Slice;
  using DimsType = BT::DimsType;
  using Box = BT::Box;
  using Num = BT::Num;
  using Coord = BT::Coord;

  constexpr virtual DimsType FindCuttingDimension(Box const& bx,
                                                  DimsType const dim) = 0;

  constexpr virtual std::pair<Box, DimsType> SwitchDimension(Slice const In,
                                                             DimsType const dim,
                                                             Box const& bx) = 0;

  constexpr virtual DimsType FindRebuildDimension(DimsType const dim) = 0;
};

template <typename Point>
struct MaxStretchDim : BaseSplitterRule<Point> {
  using BSR = BaseSplitterRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;

  void MaxStretchTag() {}

  constexpr DimsType FindCuttingDimension(
      Box const& bx, [[maybe_unused]] DimsType const dim) override {
    DimsType d(0);
    Coord diff(bx.second.pnt[0] - bx.first.pnt[0]);
    assert(Num::Geq(diff, 0));
    for (DimsType i = 1; i < BT::kDim; ++i) {
      if (Num::Gt(bx.second.pnt[i] - bx.first.pnt[i], diff)) {
        diff = bx.second.pnt[i] - bx.first.pnt[i];
        d = i;
      }
    }
    return d;
  };

  constexpr DimsType FindRebuildDimension(
      [[maybe_unused]] DimsType const dim) override {
    return 0;
  };

  constexpr std::pair<Box, DimsType> SwitchDimension(
      Slice const In, DimsType const dim,
      [[maybe_unused]] Box const& bx) override {
    return std::make_pair(BT::GetBox(In), dim);
  };
};

template <typename Point>
struct RotateDim : BaseSplitterRule<Point> {
  using BSR = BaseSplitterRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;
  using PointsIter = BT::PointsIter;

  void RotateDimTag() {}

  constexpr DimsType FindCuttingDimension([[maybe_unused]] Box const& bx,
                                          DimsType const dim) override {
    return dim;
  };

  constexpr DimsType FindRebuildDimension(
      [[maybe_unused]] DimsType const dim) override {
    return dim;
  };

  constexpr std::pair<Box, DimsType> SwitchDimension(
      Slice const In, DimsType const dim,
      [[maybe_unused]] Box const& bx) override {
    DimsType d = (dim + 1) % BT::kDim;
    // for (DimsType i = 0; i < BT::kDim; ++i, ++d) {
    //     if (!Num::Eq(In.begin()->pnt[d], std::prev(In.end())->pnt[d])) {
    //         break;
    //     }
    // }
    assert(d != dim);
    return std::make_pair(bx, d);
  };
};
}  // namespace cpdd
