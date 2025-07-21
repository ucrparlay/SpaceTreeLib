#ifndef PSPT_DEPENDENCE_SPLITTER_H
#define PSPT_DEPENDENCE_SPLITTER_H

#include <concepts>

#include "../base_tree.h"
#include "dependence/concepts.h"
#include "dependence/tree_node.h"

// NOTE: for spacial filling curve
#include "libmorton/morton.h"
#include "space_filling_curve/hilbert.c"
#include "space_filling_curve/hilbert.h"
#include "space_filling_curve/hilbert_high_dim.h"

namespace pspt {
// NOTE: ---------------- Spacial Filling Curver ---------------
template <typename Point>
struct MortonCurve {
  using BT = BaseTree<Point, MortonCurve<Point>>;
  using Slice = BT::Slice;
  using DimsType = BT::DimsType;
  using Box = BT::Box;
  using Num = BT::Num;
  using Coord = BT::Coord;
  using HyperPlane = BT::HyperPlane;

  using CurveCode = typename Point::AT::CurveCode;

  void ZTag() {}
  static std::string GetName() { return "MortonCurve"; }

  static auto Encode(Point const& p) {
    assert(std::is_integral_v<Coord>);
    // if constexpr (Point::GetDim() == 2) {
    //   return libmorton::morton2D_64_encode(p.pnt[0], p.pnt[1]);
    // } else if constexpr (Point::GetDim() == 3) {
    //   return libmorton::morton3D_64_encode(p.pnt[0], p.pnt[1], p.pnt[2]);
    // } else {
    //   static_assert("MortonCurve only supports 2D and 3D points");
    // }
    uint_fast8_t loc = 0;
    CurveCode id = 0;
    for (DimsType i = 0; i < 64 / Point::GetDim(); i++) {
      if constexpr (Point::GetDim() == 2) {
        id = id | (((p.pnt[0] >> i) & static_cast<CurveCode>(1)) << (loc++));
        id = id | (((p.pnt[1] >> i) & static_cast<CurveCode>(1)) << (loc++));
      } else if constexpr (Point::GetDim() == 3) {
        id = id | (((p.pnt[0] >> i) & static_cast<CurveCode>(1)) << (loc++));
        id = id | (((p.pnt[1] >> i) & static_cast<CurveCode>(1)) << (loc++));
        id = id | (((p.pnt[2] >> i) & static_cast<CurveCode>(1)) << (loc++));
      } else {
        for (DimsType d = 0; d < Point::GetDim(); d++) {
          id = id | (((p.pnt[d] >> i) & static_cast<CurveCode>(1)) << (loc++));
        }
      }
    }
    return id;
  }
};

template <typename Point>
struct HilbertCurve {
  using BT = BaseTree<Point, HilbertCurve<Point>>;
  using Slice = BT::Slice;
  using DimsType = BT::DimsType;
  using Box = BT::Box;
  using Num = BT::Num;
  using Coord = BT::Coord;
  using HyperPlane = BT::HyperPlane;

  using MortonCurve = MortonCurve<Point>;
  using CurveCode = typename Point::AT::CurveCode;

  void HilbertTag() {}
  static std::string GetName() { return "HilbertCurve"; }

  // TODO: optimize the encode
  static auto Encode(Point const& p) {
    assert(std::is_integral_v<Coord>);
    // auto ix = static_cast<CurveCode>(p.pnt[0]);
    // auto iy = static_cast<CurveCode>(p.pnt[1]);
    // CurveCode arr[] = {ix, iy};
    // return hilbert::hilbert_c2i(2, 32, arr);
    // return hilbert::hilbert_c2i(
    //     2, 32, reinterpret_cast<CurveCode const*>(p.GetCoords().data()));

    if constexpr (Point::GetDim() == 2) {
      return hilbert::hilbert_c2i(
          2, 32,
          reinterpret_cast<hilbert::bitmask_t const*>(p.GetCoords().data()));
      // return PSPT::LUT::hilbertXYToIndex(32, p.pnt[0], p.pnt[1]);
    } else if constexpr (Point::GetDim() == 3) {
      return PSPT::LUT::mortonToHilbert3D(MortonCurve::Encode(p), 32);
    } else {
      static_assert("HilbertCurve only supports 2D and 3D points");
    }
    // return hilbert::hilbert_c2i(2, 32, p.GetCoords().data());
  }
};

template <typename Curve>
struct SpacialFillingCurve {
  using CurveCode = typename Curve::CurveCode;

  static std::string GetSplitName() { return Curve::GetName(); }

  template <typename... Args>
  static auto Encode(Args&&... args) {
    return Curve::Encode(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto Decode(Args&&... args) {
    return Curve::Decode(std::forward<Args>(args)...);
  }
};

// NOTE: ---------------- Orthogonal Split Rule ----------------
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

  constexpr bool AllowRebuild() const { return true; };

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

  constexpr bool AllowRebuild() const { return false; };

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
struct OrthogonalSplitRule {
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
  // INFO: divide the space for binary node
  template <typename Tree, typename Slice, typename DimsType, typename Box>
  Node* DivideSpace(Tree& tree, Slice In, Slice Out, Box const& node_box,
                    Box const& input_box, DimsType dim)
    requires(IsBinaryNode<typename Tree::Interior>)
  {
    assert(Tree::WithinBox(input_box, node_box));

    // Main logic
    if (Tree::VerticalLineSplitBox(Tree::GetBoxMid(dim, node_box), input_box,
                                   dim)) {
      return tree.BuildRecursive(In, Out, dim, node_box);
    }

    auto cut_dim = dim_rule.FindCuttingDimension(node_box, dim);
    auto cut_val = Tree::GetBoxMid(cut_dim, node_box);

    // INFO: Test whether the node_box will remain unchanged after the split.
    // If the mid of box is the same as the box edge, then this time the
    // recursion will usless, the worst case is that all the mid on all
    // dimension is on the box edge, i.e., (0,1), (0,1), then a correct split
    // algorithm will handle this case
    DimsType dim_cnt = 0;
    while (dim_cnt != Tree::kDim) {
      if (!Tree::VerticalLineOnBoxEdge(cut_val, node_box, cut_dim)) {
        break;
      }
      cut_dim = dim_rule.NextDimension(cut_dim);
      cut_val = Tree::GetBoxMid(cut_dim, node_box);
      dim_cnt++;
    }

    if (dim_cnt == Tree::kDim) {  // WARN:this breaks rotation manner
      // NOTE: checks whether the node box is separatable
      return tree.BuildRecursive(In, Out, dim_rule.NextDimension(dim),
                                 node_box);
    } else if (Tree::VerticalLineSplitBox(cut_val, input_box, cut_dim)) {
      // NOTE: above while loop may changed to new dim and need to re-check
      // whether it can split the input. This is necessary as the
      // following left/right test assumes the @cut_val does not split box
      return tree.BuildRecursive(In, Out, cut_dim, node_box);
    }

    bool split_is_right = Tree::Num::Gt(cut_val, input_box.second[cut_dim]);
    assert(split_is_right || Tree::Num::Leq(cut_val, input_box.first[cut_dim]));

    typename Tree::BoxCut box_cut(
        node_box, typename Tree::Splitter(cut_val, cut_dim), split_is_right);

    Node* L = DivideSpace(tree, In, Out, box_cut.GetFirstBoxCut(), input_box,
                          dim_rule.NextDimension(cut_dim));
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

  // INFO: divide the space for multi node
  template <typename Tree, typename Slice, typename Box>
  Node* DivideSpace(Tree& tree, Slice In, Slice Out, Box const& node_box,
                    Box const& input_box)
    requires(IsMultiNode<typename Tree::Interior>)
  {
    assert(Tree::WithinBox(input_box, node_box));

    auto nodebox_split = Tree::Interior::ComputeSplitter(node_box);
    for (auto const& split : nodebox_split) {
      if (Tree::VerticalLineSplitBox(split.first, input_box, split.second)) {
        // any point on the split should be put on the right
        return tree.BuildRecursive(In, Out, node_box);
      }
    }

    if (std::ranges::all_of(nodebox_split, [&](auto const& split) {
          return Tree::VerticalLineOnBoxEdge(split.first, node_box,
                                             split.second);
        })) {  // NOTE: new box will be same as the old ones
      return tree.BuildRecursive(In, Out, node_box);
    }

    // PARA: node id contains the input box
    typename Tree::Interior::BucketType pivot_region_id = 1;
    // NOTE: considering three cases:
    // 1: the split is LESS than the left input box edge
    // 2: the split is EQ the left box edge
    // 3: same as 2 but left and right edge is same, aka. a line
    // The special judge above has prune out the case that split intersect the
    // box, therefore split Lt or Eq the RIGHT edge will cover all the cases
    for (auto const& split : nodebox_split) {
      pivot_region_id =
          pivot_region_id * 2 +
          static_cast<typename Tree::Interior::BucketType>(
              Tree::Num::Leq(split.first, input_box.second[split.second]));
    }
    pivot_region_id -= Tree::Interior::GetRegions();

    typename Tree::Interior::NodeArr tree_nodes;
    for (typename Tree::Interior::BucketType i = 0;
         i < Tree::Interior::GetRegions(); i++) {
      tree_nodes[i] = pivot_region_id == i
                          ? DivideSpace(tree, In, Out,
                                        Tree::Interior::GetBoxByRegionId(
                                            i, nodebox_split, node_box),
                                        input_box)
                          : AllocEmptyLeafNode<Slice, typename Tree::Leaf>();
    }

    return AllocInteriorNode<typename Tree::Interior>(tree_nodes, nodebox_split,
                                                      typename Tree::AugType());
  }

  // NOTE: cannot divide the points on @dim, while the points are not the same
  template <typename Tree, typename Slice, typename Box, typename... Args>
  auto HandlingUndivide(Tree& tree, Slice In, Slice Out, Box const& box,
                        Args&&... args) {
    if constexpr (IsObjectMedianSplit<PartitionRule>) {
      // NOTE: in object median, if current dimension is not divideable, then
      // switch to another dimension then continue. This works since unless all
      // points are same, otherwise we can always slice some points out.
      if constexpr (IsBinaryNode<typename Tree::Interior>) {
        // TODO: tooo brute force
        return tree.SerialBuildRecursive(
            In, Out, dim_rule.NextDimension(std::forward<Args>(args)...),
            Tree::GetBox(In));
      } else {
        return tree.BuildRecursive(In, Out, Tree::GetBox(In));
      }

    } else if constexpr (IsSpatialMedianSplit<PartitionRule>) {
      // NOTE: in spatial median, we simply reduce the box by half on
      // current dim, then switch to next dim.
      if constexpr (IsRotateDimSplit<DimRule> || IsMaxStretchDim<DimRule>) {
        return DivideSpace(tree, In, Out, box, Tree::GetBox(In),
                           std::forward<Args>(args)...);
      } else {  // define the behavior of other dim rule
        // static_assert(false);
      }
    } else {  // define the behavior of other partition rule
      // static_assert(false);
    }
  }

  DimRule const dim_rule;
  PartitionRule const partition_rule;
};
}  // namespace pspt

#endif  // PSPT_DEPENDENCE_SPLITTER_H
