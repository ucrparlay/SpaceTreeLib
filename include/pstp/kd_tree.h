#ifndef PSTP_KD_TREE_H
#define PSTP_KD_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "base_tree.h"

namespace pstp {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class KdTree
    : public BaseTree<Point, KdTree<Point, SplitRule, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio> {
 public:
  using BT = BaseTree<Point, KdTree<Point, SplitRule, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio>;

  using BucketType = typename BT::BucketType;
  using BallsType = typename BT::BallsType;
  using DimsType = typename BT::DimsType;
  using BucketSeq = typename BT::BucketSeq;
  using BallSeq = typename BT::BallSeq;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Num = Num_Comparator<Coord>;
  using Slice = typename BT::Slice;
  using Points = typename BT::Points;
  using PointsIter = typename BT::PointsIter;
  using Box = typename BT::Box;
  using BoxSeq = typename BT::BoxSeq;
  using Circle = typename BT::Circle;

  using HyperPlane = typename BT::HyperPlane;
  using HyperPlaneSeq = typename BT::HyperPlaneSeq;
  using NodeTag = typename BT::NodeTag;
  using NodeTagSeq = typename BT::NodeTagSeq;
  using NodeBox = typename BT::NodeBox;
  using NodeBoxSeq = typename BT::NodeBoxSeq;
  using Splitter = HyperPlane;
  using SplitterSeq = HyperPlaneSeq;
  using AugType = std::optional<bool>;
  using SplitRuleType = SplitRule;

  static constexpr DimsType const kMD = 2;
  using CompressNodeSplitter = std::array<HyperPlane, 1 << kMD>;

  template <uint_fast8_t kMD>
  struct KdCompressionNode;

  using Leaf =
      LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag, false>;
  struct KdInteriorNode;
  using Interior = KdInteriorNode;

  using CompressInterior = KdCompressionNode<kMD>;

  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using BoxCut = typename BT::BoxCut;

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, typename Interior, typename... Args>
  friend Node* BT::RebuildWithInsert(Node* T, Slice In, Args&&... args);

  void KdTreeTag();

  // NOTE: functions
  void Compress2Multi();

  template <typename Range>
  void Build(Range&& In);

  constexpr void DeleteTree() override;

  void BatchInsert(Slice In);

  Node* BatchInsertRecursive(Node* T, Slice In, Slice Out, DimsType d);

  // NOTE: batch delete
  // NOTE: in default, all Points to be deleted are assumed in the tree, if that
  // is not the case, using BatchDiff
  template <typename Range>
  void BatchDelete(Range&& In);

  void BatchDelete_(Slice In);

  NodeBox BatchDeleteRecursive(Node* T, Box const& bx, Slice In, Slice Out,
                               DimsType d, bool has_tomb);

  // NOTE: batch diff
  // NOTE: for the case that some Points to be deleted are not in the tree
  template <typename Range>
  void BatchDiff(Range&& In);

  void BatchDiff_(Slice In);

  // TODO: add bounding Box for batch delete recursive as well
  // WARN: fix the possible in partial deletion as well
  NodeBox BatchDiffRecursive(Node* T, Box const& bx, Slice In, Slice Out,
                             DimsType d);

  template <typename Range>
  void Flatten(Range&& Out);

  template <typename Range>
  auto KNN(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq);

  auto RangeCount(Box const& query_box);

  auto RangeCount(Circle const& cl);

  template <typename Range>
  auto RangeQuery(Box const& query_box, Range&& Out);

  void DivideRotate(Slice In, SplitterSeq& pivots, DimsType dim, BucketType idx,
                    BoxSeq& box_seq, Box const& bx);

  void PickPivots(Slice In, size_t const& n, SplitterSeq& pivots,
                  DimsType const dim, BoxSeq& box_seq, Box const& bx);

  Node* BuildRecursive(Slice In, Slice Out, DimsType dim, Box const& bx);

  Node* SerialBuildRecursive(Slice In, Slice Out, DimsType dim, Box const& bx);

  void Build_(Slice In);

  constexpr static char const* GetTreeName() { return "KdTree"; }

  SplitRule split_rule_;
};

}  // namespace pstp

#include "kd_tree_impl/kd_batch_delete.hpp"
#include "kd_tree_impl/kd_batch_diff.hpp"
#include "kd_tree_impl/kd_batch_insert.hpp"
#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"

#endif  // PSTP_KD_TREE_H
