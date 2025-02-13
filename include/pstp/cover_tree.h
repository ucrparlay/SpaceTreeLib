#ifndef PSTP_COVER_TREE_H
#define PSTP_COVER_TREE_H

#include <functional>
#include <utility>

#include "base_tree.h"

namespace pstp {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class CoverTree
    : public BaseTree<Point, CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio> {
 public:
  using BT = BaseTree<Point, CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio>;

  using BucketType = BT::BucketType;
  using BallsType = BT::BallsType;
  using BucketSeq = BT::BucketSeq;
  using DepthType = BT::DepthType;
  using BallSeq = BT::BallSeq;
  using DimsType = BT::DimsType;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Num = Num_Comparator<Coord>;
  using Slice = BT::Slice;
  using Points = BT::Points;
  using PointsIter = BT::PointsIter;
  using Box = BT::Box;
  using BoxSeq = BT::BoxSeq;
  using Circle = BT::Circle;

  using HyperPlane = BT::HyperPlane;
  using HyperPlaneSeq = BT::HyperPlaneSeq;
  using Splitter =
      parlay::sequence<Point>;  // every center of circle in one
                                // level marks a region, the radius can be
                                // passed during the traversal
  using SplitterSeq = parlay::sequence<Splitter>;
  using NodeTagSeq = BT::NodeTagSeq;
  using NodeBoxSeq = BT::NodeBoxSeq;
  using NodeBoolean = BT::NodeBoolean;
  using NodeBox = BT::NodeBox;
  using AugType = std::optional<bool>;

  struct CoverInteriorNode;

  using SplitRuleType = SplitRule;
  using Leaf =
      LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag, false>;
  using Interior = CoverInteriorNode;
  using CoverNodeArr = Interior::CoverNodeArr;
  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using BoxCut = typename BT::BoxCut;
  using CoverCircle = BT::CoverCircle;

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, typename Interior, typename... Args>
  friend Node* BT::RebuildWithInsert(Node* T, Slice In, Args&&... args);

  void CoverTreeTag();

  // NOTE: functions
  template <typename Range, typename... Args>
  void Build(Range&& In, Args&&... args);

  template <typename Range>
  void BatchInsert(Range&& In);

  NodeBoolean PointInsertRecursive(Node* T, Point const& center, Point const& p,
                                   DepthType deep);

  void BuildUpwards(CoverCircle& root_cc, Point const& p);

  constexpr void DeleteTree() override;

  // NOTE: batch delete
  // NOTE: in default, all Points to be deleted are assumed in the tree
  template <typename Range>
  void BatchDelete(Range&& In);

  void BatchDelete_(Slice In);

  Node* BatchDeleteRecursive(Node* T, Slice In, Slice Out, Box const& box,
                             bool has_tomb);

  // NOTE: batch diff
  // NOTE: for the case that some Points to be deleted are not in the tree
  template <typename Range>
  void BatchDiff(Range&& In);

  template <typename Range>
  void Flatten(Range&& Out);

  template <typename Range>
  auto KNN(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq);

  auto RangeCount(Box const& query_box);

  auto RangeCount(Circle const& cl);

  template <typename Range>
  auto RangeQuery(Box const& query_box, Range&& Out);

  void Build_(Slice In);

  void Build_(Slice In, Box const& box);

  void SerialSplit(Slice In, DimsType dim, DimsType idx, Box const& box,
                   parlay::sequence<BallsType>& sums);

  void SerialSplitSkeleton(Node* T, Slice In, DimsType dim, DimsType idx,
                           parlay::sequence<BallsType>& sums);

  void DivideRotate(HyperPlaneSeq& pivots, DimsType dim, BucketType idx,
                    BoxSeq& box_seq, Box const& box);

  void PickPivots(Slice In, size_t const& n, HyperPlaneSeq& pivots,
                  DimsType const dim, BoxSeq& box_seq, Box const& bx);

  Node* BuildRecursive(Slice In, Slice Out, Box const& bx);

  Node* SerialBuildRecursive(Slice In, Slice Out, Box const& bx,
                             bool checked_duplicate);

  void BatchInsert_(Slice In);

  Node* BatchInsertRecursive(Node* T, Slice In, Slice Out);

  void BatchDiff_(Slice In);

  Node* BatchDiffRecursive(Node* T, Slice In, Slice Out);

  constexpr static char const* GetTreeName() { return "CoverTree"; }

  SplitRule split_rule_;
  CoverCircle root_cover_circle_;
  size_t alloc_dummy_num_ = 0;
  size_t alloc_empty_num_ = 0;
  size_t alloc_normal_num_ = 0;
  size_t alloc_interior_num_ = 0;
};

}  // namespace pstp

// #include "cover_tree_impl/cover_batch_delete.hpp"
// #include "cover_tree_impl/cover_batch_diff.hpp"
#include "cover_tree_impl/cover_batch_insert.hpp"
#include "cover_tree_impl/cover_build_tree.hpp"
#include "cover_tree_impl/cover_inter_node.hpp"
#include "cover_tree_impl/cover_override.hpp"

#endif  // PSTP_COVER_TREE_H
