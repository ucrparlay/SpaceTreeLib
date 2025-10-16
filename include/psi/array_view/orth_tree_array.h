#ifndef PSI_ARRAY_VIEW_ORTH_TREE_ARRAY_H_
#define PSI_ARRAY_VIEW_ORTH_TREE_ARRAY_H_

#include <array>
#include <functional>
#include <utility>
#include <vector>

#include "base_tree_array.h"

namespace psi {
namespace array_view {

//==============================================================================
// ARRAY-BASED ORTHOGONAL TREE
//
// Array-based implementation of Orthogonal Tree (multi-way tree).
// Provides the same interface as pointer-based OrthTree with improved
// cache locality and reduced memory overhead.
//==============================================================================

template <typename Point, typename SplitRule, typename LeafAugType,
          typename InteriorAugType, uint_fast8_t kMD = 2,
          uint_fast8_t kSkHeight = 6, uint_fast8_t kImbaRatio = 30>
class OrthTreeArray
    : public BaseTreeArray<
          Point,
          OrthTreeArray<Point, SplitRule, LeafAugType, InteriorAugType, kMD,
                        kSkHeight, kImbaRatio>,
          kSkHeight, kImbaRatio> {
 public:
  static constexpr size_t kSplitterNum = kMD;
  static constexpr size_t kNodeRegions = 1 << kMD;

  using BT =
      BaseTreeArray<Point,
                    OrthTreeArray<Point, SplitRule, LeafAugType,
                                  InteriorAugType, kMD, kSkHeight, kImbaRatio>,
                    kSkHeight, kImbaRatio>;

  using NodeIndex = typename BT::NodeIndex;
  using BucketType = typename BT::BucketType;
  using BallsType = typename BT::BallsType;
  using BucketSeq = typename BT::BucketSeq;
  using BallSeq = typename BT::BallSeq;
  using DimsType = typename BT::DimsType;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Num = Num_Comparator<Coord>;
  using Slice = typename BT::Slice;
  using Points = typename BT::Points;
  using PointsIter = typename BT::PointsIter;
  using Box = typename BT::Box;
  using BoxSeq = typename BT::BoxSeq;

  using HyperPlane = typename BT::HyperPlane;
  using HyperPlaneSeq = typename BT::HyperPlaneSeq;
  using Splitter = std::array<HyperPlane, kSplitterNum>;
  using SplitterSeq = parlay::sequence<Splitter>;
  using SplitRuleType = SplitRule;

  //============================================================================
  // COMPACT MULTI-WAY NODE
  //============================================================================

  struct CompactMultiNode {
    std::array<NodeIndex, kNodeRegions> children;  // Child indices
    Splitter split;                                // Splitting hyperplanes
    InteriorAugType aug;                           // Augmentation
    uint32_t size;                                 // Subtree size
    uint8_t is_leaf;                               // Leaf flag

    // For leaf nodes
    uint32_t leaf_start_idx;
    uint32_t leaf_count;

    CompactMultiNode()
        : size(0), is_leaf(false), leaf_start_idx(0), leaf_count(0) {
      children.fill(BT::NULL_INDEX);
    }
  };

  using Node = CompactMultiNode;

  //============================================================================
  // STORAGE
  //============================================================================

  std::vector<CompactMultiNode> nodes_;
  std::vector<Point> leaf_points_;

  //============================================================================
  // PUBLIC API
  //============================================================================

  void OrthTreeArrayTag() {}

  template <typename Range, typename... Args>
  void Build(Range&& In, Args&&... args);

  template <typename Range>
  void BatchInsert(Range&& In);

  constexpr void DeleteTree() override;

  template <typename Range>
  void BatchDelete(Range&& In);

  template <typename Range>
  void BatchDiff(Range&& In);

  template <typename Range>
  void Flatten(Range&& Out);

  template <typename Range>
  auto KNN(NodeIndex idx, Point const& q, kBoundedQueue<Point, Range>& bq);

  auto RangeCount(Box const& query_box);

  template <typename Range>
  auto RangeQuery(Box const& query_box, Range&& Out);

  constexpr static char const* GetTreeName() { return "OrthTreeArray"; }

  constexpr static char const* CheckHasBox() {
    if constexpr (HasBox<InteriorAugType>)
      return "HasBox";
    else
      return "NoBox";
  }

 private:
  void Build_(Slice In);
  void Build_(Slice In, Box const& box);

  NodeIndex BuildRecursive(Slice In, Slice Out, Box const& box);

  NodeIndex SerialBuildRecursive(Slice In, Slice Out, Box const& box,
                                 bool checked_duplicate);

  // Node allocation
  NodeIndex AllocNode();
  NodeIndex AllocLeafNode(Slice In);
  void FreeNode(NodeIndex idx);

  // Node access
  CompactMultiNode& GetNode(NodeIndex idx) { return nodes_[idx]; }
  CompactMultiNode const& GetNode(NodeIndex idx) const { return nodes_[idx]; }

  bool IsLeaf(NodeIndex idx) const {
    return idx != BT::NULL_INDEX && nodes_[idx].is_leaf;
  }

  SplitRule split_rule_;
  bool fixed_box = false;
};

}  // namespace array_view
}  // namespace psi

#include "orth_tree_array_impl/orth_batch_delete.hpp"
#include "orth_tree_array_impl/orth_batch_diff.hpp"
#include "orth_tree_array_impl/orth_batch_insert.hpp"
#include "orth_tree_array_impl/orth_build_tree.hpp"
#include "orth_tree_array_impl/orth_override.hpp"

#endif  // PSI_ARRAY_VIEW_ORTH_TREE_ARRAY_H_
