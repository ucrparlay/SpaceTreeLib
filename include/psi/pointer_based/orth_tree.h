#ifndef PSI_POINTER_BASED_ORTH_TREE_H_
#define PSI_POINTER_BASED_ORTH_TREE_H_

#include <functional>
#include <utility>

#include "base_tree.h"

namespace psi {

template <typename TypeTrait>
class OrthTree : public BaseTree<TypeTrait, OrthTree<TypeTrait>> {
 public:
  static constexpr uint_fast8_t kMD = TypeTrait::Point::GetDim();
  static constexpr size_t kSplitterNum = kMD;
  static constexpr size_t kNodeRegions = 1 << kMD;

  using BT = BaseTree<TypeTrait, OrthTree<TypeTrait>>;
  using Geo = GeoBase<TypeTrait>;

  using Point = typename BT::Point;
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
  // using Circle = BT::Circle;

  using HyperPlane = BT::HyperPlane;
  using HyperPlaneSeq = BT::HyperPlaneSeq;
  using Splitter = std::array<HyperPlane, kSplitterNum>;
  using SplitterSeq = parlay::sequence<Splitter>;
  using NodeTagSeq = BT::NodeTagSeq;
  using NodeBoxSeq = BT::NodeBoxSeq;
  using NodeBox = BT::NodeBox;
  // using AugType = std::optional<bool>;

  using LeafAugType = typename TypeTrait::LeafAugType;
  using InteriorAugType = typename TypeTrait::InteriorAugType;

  struct OrthInteriorNode;

  using SplitRule = typename TypeTrait::SplitRule;
  using SplitRuleType = typename TypeTrait::SplitRule;
  using Leaf = LeafNode<Point, Slice, BT::kLeaveWrap, LeafAugType,
                        parlay::move_assign_tag>;
  using Interior = OrthInteriorNode;
  using OrthNodeArr = Interior::OrthNodeArr;
  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using BoxCut = typename BT::BoxCut;

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, typename Interior, typename PrepareFunc,
            typename... Args>
  friend Node* BT::RebuildWithInsert(Node* T, PrepareFunc prepare_func,
                                     Slice In, Args&&... args);

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

  // auto RangeCount(Circle const& cl);

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

  Node* BatchInsertRecursive(Node* T, Slice In, Slice Out, Box const& bx);

  void BatchDiff_(Slice In);

  void SetBoundingBox(Box const& box) {
    this->tree_box_ = box;
    fixed_box = true;
  }

  Node* BatchDiffRecursive(Node* T, Slice In, Slice Out);

  constexpr static char const* GetTreeName() { return "OrthTree"; }
  constexpr static char const* CheckHasBox() {
    if constexpr (HasBox<InteriorAugType>)
      return "HasBox";
    else
      return "NoBox";
  }

  SplitRule split_rule_;
  bool fixed_box = false;
  size_t alloc_dummy_num_ = 0;
  size_t alloc_empty_num_ = 0;
  size_t alloc_normal_num_ = 0;
  size_t alloc_interior_num_ = 0;
};

}  // namespace psi

#include "orth_tree_impl/orth_batch_delete.hpp"
#include "orth_tree_impl/orth_batch_diff.hpp"
#include "orth_tree_impl/orth_batch_insert.hpp"
#include "orth_tree_impl/orth_build_tree.hpp"
#include "orth_tree_impl/orth_inter_node.hpp"
#include "orth_tree_impl/orth_override.hpp"

#endif  // PSI_POINTER_BASED_ORTH_TREE_H_