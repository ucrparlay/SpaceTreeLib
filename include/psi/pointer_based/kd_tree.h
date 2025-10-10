#ifndef PSI_POINTER_BASED_KD_TREE_H_
#define PSI_POINTER_BASED_KD_TREE_H_

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "../dependence/concepts.h"
#include "base_tree.h"

namespace psi {
namespace pointer_based {

template <typename TypeTrait>
class KdTree : public BaseTree<TypeTrait, KdTree<TypeTrait>> {
 public:
  using BT = BaseTree<TypeTrait, KdTree<TypeTrait>>;
  using Geo = GeoBase<TypeTrait>;

  using Point = typename TypeTrait::Point;
  using BucketType = typename TypeTrait::BucketType;
  using BallsType = typename TypeTrait::BallsType;
  using DimsType = typename TypeTrait::DimsType;
  using BucketSeq = typename TypeTrait::BucketSeq;
  using BallSeq = typename TypeTrait::BallSeq;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Num = Num_Comparator<Coord>;
  using Slice = typename TypeTrait::Slice;
  using Points = typename TypeTrait::Points;
  using PointsIter = typename TypeTrait::PointsIter;
  using Box = typename TypeTrait::Box;
  using BoxSeq = typename TypeTrait::BoxSeq;

  using HyperPlane = typename TypeTrait::HyperPlane;
  using HyperPlaneSeq = typename TypeTrait::HyperPlaneSeq;
  using NodeTag = typename TypeTrait::NodeTag;
  using NodeTagSeq = typename TypeTrait::NodeTagSeq;
  using NodeBox = typename TypeTrait::NodeBox;
  using NodeBoxSeq = typename TypeTrait::NodeBoxSeq;
  using Splitter = HyperPlane;
  using SplitterSeq = HyperPlaneSeq;
  using Circle = BT::NormalCircle;

  using SplitRuleType = TypeTrait::SplitRule;
  using SplitRule = typename TypeTrait::SplitRule;
  using LeafAugType = typename TypeTrait::LeafAugType;
  using InteriorAugType = typename TypeTrait::InteriorAugType;

  static constexpr DimsType const kMD = 2;
  using CompressNodeSplitter = std::array<HyperPlane, 1 << kMD>;

  template <uint_fast8_t kMD>
  struct KdCompressionNode;

  using Leaf = LeafNode<Point, Slice, BT::kLeaveWrap, LeafAugType,
                        parlay::move_assign_tag>;
  struct KdInteriorNode;
  using Interior = KdInteriorNode;

  using CompressInterior = KdCompressionNode<kMD>;

  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using BoxCut = typename BT::BoxCut;

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, typename Interior, typename PrepareFunc,
            typename... Args>
  friend Node* BT::RebuildWithInsert(Node* T, PrepareFunc prepare_func,
                                     Slice In, Args&&... args);

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

  constexpr static char const* CheckHasBox() {
    if constexpr (HasBox<InteriorAugType>)
      return "HasBox";
    else
      return "NoBox";
  }

  SplitRule split_rule_;
};

}  // namespace pointer_based
}  // namespace psi

#include "kd_tree_impl/kd_batch_delete.hpp"
#include "kd_tree_impl/kd_batch_diff.hpp"
#include "kd_tree_impl/kd_batch_insert.hpp"
#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"

#endif  // PSI_POINTER_BASED_KD_TREE_H_