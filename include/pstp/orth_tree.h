#ifndef PSTP_ORTH_TREE_H
#define PSTP_ORTH_TREE_H

#include <functional>
#include <utility>

#include "base_tree.h"

namespace pstp {

template <typename Point, typename SplitRule, uint_fast8_t kMD = 2,
          uint_fast8_t kSkHeight = 6, uint_fast8_t kImbaRatio = 30>
class OrthTree
    : public BaseTree<Point,
                      OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio> {
 public:
  static constexpr size_t kSplitterNum = kMD;
  static constexpr size_t kNodeRegions = 1 << kMD;

  using BT =
      BaseTree<Point, OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>,
               kSkHeight, kImbaRatio>;

  using BucketType = BT::BucketType;
  using BallsType = BT::BallsType;
  using BucketSeq = BT::BucketSeq;
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
  using Splitter = std::array<HyperPlane, kSplitterNum>;
  using SplitterSeq = parlay::sequence<Splitter>;
  using NodeTagSeq = BT::NodeTagSeq;
  using NodeBoxSeq = BT::NodeBoxSeq;
  using NodeBox = BT::NodeBox;
  using AugType = std::optional<bool>;

  // struct KdInteriorNode;
  struct OrthInteriorNode;

  using SplitRuleType = SplitRule;
  using Leaf =
      LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag, false>;
  using Interior = OrthInteriorNode;
  using OrthNodeArr = Interior::OrthNodeArr;
  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using BoxCut = typename BT::BoxCut;

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, typename Interior, typename... Args>
  friend Node* BT::RebuildWithInsert(Node* T, Slice In, Args&&... args);

  void OrthTreeTag();

  // NOTE: functions
  template <typename Range, typename... Args>
  void Build(Range&& In, Args&&... args);

  template <typename Range>
  void BatchInsert(Range&& In);

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

  constexpr static char const* GetTreeName() { return "OrthTree"; }

  SplitRule split_rule_;
  size_t alloc_dummy_num_ = 0;
  size_t alloc_empty_num_ = 0;
  size_t alloc_normal_num_ = 0;
  size_t alloc_interior_num_ = 0;
};

}  // namespace pstp

#include "orth_tree_impl/orth_batch_delete.hpp"
#include "orth_tree_impl/orth_batch_diff.hpp"
#include "orth_tree_impl/orth_batch_insert.hpp"
#include "orth_tree_impl/orth_build_tree.hpp"
#include "orth_tree_impl/orth_inter_node.hpp"
#include "orth_tree_impl/orth_override.hpp"

#endif  // PSTP_ORTH_TREE_H
