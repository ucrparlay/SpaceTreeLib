#ifndef PSTP_DEPENDENCE_SPLITTER_H
#define PSTP_DEPENDENCE_SPLITTER_H

#include "../base_tree.h"
#include "dependence/concepts.h"
#include "dependence/tree_node.h"

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

  constexpr virtual DimsType const FindRebuildDimension(
      DimsType const dim) const = 0;

  constexpr virtual DimsType const NextDimension(DimsType const dim) const = 0;
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
  static std::string GetName() { return "MaxStretchDim"; }

  constexpr DimsType const FindCuttingDimension(Box const& bx,
                                                DimsType) const override {
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

  constexpr DimsType const FindRebuildDimension(DimsType) const override {
    return 0;
  };

  // TODO: this is wired as the next dimension should return the dimension with
  // second largest
  constexpr DimsType const NextDimension(DimsType) const override { return 0; };
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
  static std::string GetName() { return "RotateDim"; }

  constexpr DimsType const FindCuttingDimension(
      [[maybe_unused]] Box const& bx, DimsType const dim) const override {
    return dim;
  };

  constexpr DimsType const FindRebuildDimension(
      [[maybe_unused]] DimsType const dim) const override {
    return dim;
  };

  constexpr DimsType const NextDimension(DimsType dim) const override {
    return (dim + 1) % BT::kDim;
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
  using IterHyperPair = std::pair<PointsIter, std::optional<HyperPlane>>;

  constexpr virtual IterHyperPair const SplitInput(Slice In, DimsType const dim,
                                                   Box const& box) const = 0;
  constexpr virtual HyperPlane const SplitSample(Slice In, DimsType const dim,
                                                 Box const& box) const = 0;

  constexpr virtual bool AllowRebuild() const = 0;
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
  using IterHyperPair = BSR::IterHyperPair;

  void ObjectMedianTag() {}

  static std::string GetName() { return "ObjectMedian"; }

  constexpr bool AllowRebuild() { return true; };

  constexpr IterHyperPair const SplitInput(
      Slice In, DimsType const dim,
      [[maybe_unused]] Box const& box) const override {
    size_t n = In.size();
    std::ranges::nth_element(In.begin(), In.begin() + n / 2, In.end(),
                             [&](Point const& p1, Point const& p2) {
                               return Num::Lt(p1.pnt[dim], p2.pnt[dim]);
                             });

    auto split_iter =
        std::ranges::partition(In.begin(), In.begin() + n / 2,
                               [&](Point const& p) {
                                 return Num::Lt(p.pnt[dim], In[n / 2].pnt[dim]);
                               })
            .begin();

    if (split_iter == In.begin()) {  // NOTE: handle duplicated medians
      split_iter =
          std::ranges::partition(In.begin() + n / 2, In.end(),
                                 [&](Point const& p) {
                                   return Num::Eq(p.pnt[dim],
                                                  In[n / 2].pnt[dim]);
                                 })
              .begin();  // NOTE: now all duplicated median is on the left
    }
    // return split_iter;
    if (split_iter <= In.begin() + n / 2) {  // NOTE: split is on left half
      return IterHyperPair(split_iter, HyperPlane(In[n / 2].pnt[dim], dim));
    } else if (split_iter != In.end()) {  // NOTE: split is on right half
      auto min_elem_iter = std::ranges::min_element(
          split_iter, In.end(), [&](Point const& p1, Point const& p2) {
            return Num::Lt(p1.pnt[dim], p2.pnt[dim]);
          });
      return IterHyperPair(split_iter,
                           HyperPlane(min_elem_iter->pnt[dim], dim));
    } else {  // NOTE: all the same
      return IterHyperPair(split_iter, std::nullopt);
    }
  }

  // TODO: the sample may handle the duplicates as well
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
  using IterHyperPair = BSR::IterHyperPair;

  void SpatialMedianTag() {}

  static std::string GetName() { return "SpatialMedian"; }

  constexpr bool AllowRebuild() { return false; };

  constexpr IterHyperPair const SplitInput(Slice In, DimsType const dim,
                                           Box const& box) const override {
    auto split_iter = std::ranges::partition(In, [&](Point const& p) {
                        return Num::Lt(p.pnt[dim], BT::GetBoxMid(dim, box));
                      }).begin();
    if (split_iter == In.begin() || split_iter == In.end()) {
      return IterHyperPair(split_iter, std::nullopt);
    } else {
      return IterHyperPair(split_iter,
                           HyperPlane(BT::GetBoxMid(dim, box), dim));
    }
  }

  constexpr HyperPlane const SplitSample([[maybe_unused]] Slice In,
                                         DimsType const dim,
                                         Box const& box) const override {
    return HyperPlane(BT::GetBoxMid(dim, box), dim);
  }
};

template <class DimRule, class PartitionRule>
struct SplitRule {
  using DimRuleType = DimRule;
  using PartitionRuleType = PartitionRule;

  static std::string GetSplitName() {
    return DimRule::GetName() + "-" + PartitionRule::GetName();
  }

  // NOTE: dimension
  template <typename... Args>
  auto FindCuttingDimension(Args&&... args) {
    return dim_rule.FindCuttingDimension(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto FindRebuildDimension(Args&&... args) {
    return dim_rule.FindRebuildDimension(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto NextDimension(Args&&... args) {
    return dim_rule.NextDimension(std::forward<Args>(args)...);
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

  // NOTE: query whether to launch the rebuild
  template <typename... Args>
  auto AllowRebuild(Args&&... args) {
    return partition_rule.AllowRebuild(std::forward<Args>(args)...);
  }

  // NOTE: helper for handling the duplicate
  // NOTE: divide the space until the split cut the input box
  template <typename Tree, typename Slice, typename DimsType, typename Box>
  Node* DivideSpace(Tree& tree, Slice In, Slice Out, DimsType dim,
                    Box const& node_box, Box const& input_box) {
    assert(Tree::WithinBox(input_box, node_box));

    // the mid of the box split the box for input points
    if (Tree::VerticalLineIntersectBoxExclude(Tree::GetBoxMid(dim, node_box),
                                              input_box, dim)) {
      return tree.BuildRecursive(In, Out, dim, node_box);
    }

    auto cut_dim = dim_rule.FindCuttingDimension(node_box, dim);
    auto cut_val = Tree::GetBoxMid(cut_dim, node_box);

    // if the line is on the left of the boundary, then we should generate the
    // box to the right, and vice versa.
    // otherwise determine the position by comparing the mid of the box
    // This should work for the max stretch dim as well
    bool split_is_right = Tree::Num::Geq(Tree::GetBoxMid(cut_dim, node_box),
                                         input_box.second[cut_dim]);
    assert(split_is_right || Tree::Num::Leq(Tree::GetBoxMid(cut_dim, node_box),
                                            input_box.first[cut_dim]));

    typename Tree::BoxCut box_cut(
        node_box, typename Tree::Splitter(cut_val, cut_dim), split_is_right);

    Node* L = DivideSpace(tree, In, Out, dim_rule.NextDimension(cut_dim),
                          box_cut.GetFirstBoxCut(), input_box);
    Node* R = AllocEmptyLeafNode<Slice, typename Tree::Leaf>();
    assert(Tree::WithinBox(input_box, box_cut.GetBox()));

    if (!split_is_right) {
      assert(Tree::Num::Leq(Tree::GetBoxMid(cut_dim, node_box),
                            input_box.first[cut_dim]));
      std::ranges::swap(L, R);
    }

    return AllocInteriorNode<typename Tree::Interior>(
        L, R, box_cut.GetHyperPlane(), typename Tree::AugType());
  }

  // NOTE: cannot divide the points on @dim, while the points are not the same
  template <typename Tree, typename Slice, typename DimsType, typename Box>
  auto HandlingUndivide(Tree& tree, Slice In, Slice Out, DimsType dim,
                        Box const& box) {
    if constexpr (IsObjectMedianSplit<PartitionRule>) {
      // NOTE: in object median, if current dimension is not divideable, then
      // switch to another dimension then continue. This works since unless all
      // points are same, otherwise we can always slice some points out.
      return tree.SerialBuildRecursive(In, Out, dim_rule.NextDimension(dim),
                                       Tree::GetBox(In));

    } else if constexpr (IsSpatialMedianSplit<PartitionRule>) {
      // NOTE: in spatial median, we simply reduce the box by half on
      // current dim, then switch to next dim.
      if constexpr (IsRotateDimSplit<DimRule> || IsMaxStretchDim<DimRule>) {
        // TODO: add support for multi node
        assert(IsBinaryNode<typename Tree::Interior>);
        return DivideSpace(tree, In, Out, dim, box, Tree::GetBox(In));
      } else {  // define the behavior of other dim rule
        static_assert(false);
      }
    } else {  // define the behavior of other partition rule
      static_assert(false);
    }
  }

  DimRule const dim_rule;
  PartitionRule const partition_rule;
};

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_SPLITTER_H
