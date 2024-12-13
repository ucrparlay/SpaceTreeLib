#ifndef PSTP_DEPENDENCE_SPLITTER_H
#define PSTP_DEPENDENCE_SPLITTER_H

#include "../base_tree.h"

namespace pstp {
template <typename Point>
struct BaseSplitDimRule {
  using BT = BaseTree<Point, BaseSplitDimRule<Point>>;
  using Slice = BT::Slice;
  using DimsType = BT::DimsType;
  using Box = BT::Box;
  using Num = BT::Num;
  using Coord = BT::Coord;
  using HyperPlane = BT::HyperPlane;

  constexpr virtual DimsType const FindCuttingDimension(
      Box const& bx, DimsType const dim) const = 0;

  constexpr virtual std::pair<Box, DimsType> const SwitchDimension(
      Slice const In, DimsType const dim, Box const& bx) const = 0;

  constexpr virtual DimsType const FindRebuildDimension(
      DimsType const dim) const = 0;

  // TODO: the spliiter should deterine how to split as well
};

template <typename Point>
struct MaxStretchDim : BaseSplitDimRule<Point> {
  using BSR = BaseSplitDimRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;
  using HyperPlane = BT::HyperPlane;

  void MaxStretchTag() {}

  constexpr DimsType const FindCuttingDimension(
      Box const& bx, [[maybe_unused]] DimsType const dim) const override {
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

  constexpr DimsType const FindRebuildDimension(
      [[maybe_unused]] DimsType const dim) const override {
    return 0;
  };

  constexpr std::pair<Box, DimsType> const SwitchDimension(
      Slice const In, DimsType const dim,
      [[maybe_unused]] Box const& bx) const override {
    return std::make_pair(BT::GetBox(In), dim);
  };
};

template <typename Point>
struct RotateDim : BaseSplitDimRule<Point> {
  using BSR = BaseSplitDimRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;
  using PointsIter = BT::PointsIter;
  using HyperPlane = BT::HyperPlane;

  void RotateDimTag() {}

  constexpr DimsType const FindCuttingDimension(
      [[maybe_unused]] Box const& bx, DimsType const dim) const override {
    return dim;
  };

  constexpr DimsType const FindRebuildDimension(
      [[maybe_unused]] DimsType const dim) const override {
    return dim;
  };

  constexpr std::pair<Box, DimsType> const SwitchDimension(
      [[maybe_unused]] Slice const In, DimsType const dim,
      [[maybe_unused]] Box const& bx) const override {
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

template <typename Point>
struct BaseSplitPartitionRule {
  using BT = BaseTree<Point, BaseSplitPartitionRule<Point>>;
  using Slice = BT::Slice;
  using DimsType = BT::DimsType;
  using Box = BT::Box;
  using Num = BT::Num;
  using Coord = BT::Coord;
  using HyperPlane = BT::HyperPlane;
  using PointsIter = BT::PointsIter;

  constexpr virtual PointsIter const SplitInput(Slice In, DimsType const dim,
                                                Box const& box) const = 0;
  constexpr virtual HyperPlane const SplitSample(Slice In, DimsType const dim,
                                                 Box const& box) const = 0;
};

template <typename Point>
struct ObjectMedian : BaseSplitPartitionRule<Point> {
  using BSR = BaseSplitPartitionRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;
  using PointsIter = BSR::PointsIter;
  using HyperPlane = BSR::HyperPlane;

  void ObjectMedianTag() {}

  constexpr PointsIter const SplitInput(
      Slice In, DimsType const dim,
      [[maybe_unused]] Box const& box) const override {
    size_t n = In.size();
    std::ranges::nth_element(In.begin(), In.begin() + n / 2, In.end(),
                             [&](Point const& p1, Point const& p2) {
                               return Num::Lt(p1.pnt[dim], p2.pnt[dim]);
                             });

    std::ranges::subrange _2ndGroup = std::ranges::partition(
        In.begin(), In.begin() + n / 2, [&](Point const& p) {
          return Num::Lt(p.pnt[dim], In[n / 2].pnt[dim]);
        });

    if (_2ndGroup.begin() == In.begin()) {  // NOTE: handle duplicated medians
      _2ndGroup = std::ranges::partition(
          In.begin() + n / 2, In.end(), [&](Point const& p) {
            return Num::Eq(p.pnt[dim], In[n / 2].pnt[dim]);
          });  // NOTE: now all duplicated median is on the left
    }
    return _2ndGroup.begin();
  }

  constexpr HyperPlane const SplitSample(
      Slice In, DimsType const dim,
      [[maybe_unused]] Box const& box) const override {
    size_t n = In.size();
    std::ranges::nth_element(In, In.begin() + n / 2,
                             [&](Point const& p1, Point const& p2) {
                               return Num::Lt(p1.pnt[dim], p2.pnt[dim]);
                             });
    return HyperPlane(In[n / 2][dim], dim);
  }
};

template <typename Point>
struct SpatialMedian : BaseSplitPartitionRule<Point> {
  using BSR = BaseSplitPartitionRule<Point>;
  using BT = BSR::BT;
  using Slice = BSR::Slice;
  using DimsType = BSR::DimsType;
  using Box = BSR::Box;
  using Coord = BSR::Coord;
  using Num = BSR::Num;
  using PointsIter = BSR::PointsIter;
  using HyperPlane = BSR::HyperPlane;

  void SpatialMedianTag() {}

  constexpr PointsIter const SplitInput(Slice In, DimsType const dim,
                                        Box const& box) const override {
    return std::ranges::partition(In,
                                  [&](Point const& p) {
                                    return Num::Lt(p.pnt[dim],
                                                   BT::GetBoxMid(dim, box));
                                  })
        .begin();
  }

  constexpr HyperPlane const SplitSample([[maybe_unused]] Slice In,
                                         DimsType const dim,
                                         Box const& box) const override {
    return HyperPlane(BT::GetBoxMid(dim, box), dim);
  }
};

template <class DimRule, class PartitionRule>
struct SplitRule {
  // NOTE: dimension
  template <typename... Args>
  auto FindCuttingDimension(Args&&... args) {
    return dim_rule.FindCuttingDimension(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto SwitchDimension(Args&&... args) {
    return dim_rule.SwitchDimension(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto FindRebuildDimension(Args&&... args) {
    return dim_rule.FindRebuildDimension(std::forward<Args>(args)...);
  }

  // NOTE: serial parititon used in algorithm
  template <typename... Args>
  auto SplitInput(Args&&... args) {
    return partition_rule.SplitInput(std::forward<Args>(args)...);
  }

  // NOTE: split the sample in order to get the hyperplane
  template <typename... Args>
  auto SplitSample(Args&&... args) {
    return partition_rule.SplitSample(std::forward<Args>(args)...);
  }

  DimRule const dim_rule;
  PartitionRule const partition_rule;
};

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_SPLITTER_H
