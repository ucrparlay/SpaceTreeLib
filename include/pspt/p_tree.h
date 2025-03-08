#ifndef PSPT_P_TREE_H
#define PSPT_P_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "base_tree.h"

namespace pspt {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class PTree
    : public BaseTree<Point, PTree<Point, SplitRule, kSkHeight, kImbaRatio>,
                      kSkHeight, kImbaRatio> {
 public:
  using BT = BaseTree<Point, PTree<Point, SplitRule, kSkHeight, kImbaRatio>,
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
  using Circle = typename BT::NormalCircle;

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

  // NOTE: for the CPAM
  using EncodedType = typename BT::EncodedType;
  using IDType = typename BT::IDType;
  using CpamKey = pair<EncodedType, IDType>;  // morton_id, id
  using CpamVal = Point;
  using CpamAug = pair<Box, size_t>;

  struct CpamEntry {
    using key_t = CpamKey;
    using val_t = CpamVal;
    using aug_t = CpamAug;

    static inline bool comp(key_t a, key_t b) { return a < b; }
    static aug_t get_empty() { return make_pair(BT::GetEmptyBox(), 0); }
    static aug_t from_entry([[maybe_unused]] key_t k, val_t v) {
      return make_pair(Box(v, v), 1);
    }
    static aug_t combine(aug_t a, aug_t b) {
      return make_pair(BT::GetBox(a.first, b.first), a.second + b.second);
    }
  };

  using CpamAugMap = cpam::aug_map<CpamEntry, 32>;
  using par = std::tuple<typename CpamEntry::key_t, typename CpamEntry::val_t>;

  using Leaf = CpamAugMap::Tree::node;
  using Interior = CpamAugMap::Tree::node;

  // NOTE: general tree structure
  void PTreeTag();

  template <typename Range>
  void Build(Range&& In);

  void Build_(Slice In);

  constexpr void DeleteTree() override;

  void BatchInsert(Slice In);

  void BatchInsert_(Slice In);

  // NOTE: batch delete
  // NOTE: in default, all Points to be deleted are assumed in the tree, if that
  // is not the case, using BatchDiff
  template <typename Range>
  void BatchDelete(Range&& In);

  void BatchDelete_(Slice In);

  // NOTE: batch diff
  // NOTE: for the case that some Points to be deleted are not in the tree
  template <typename Range>
  void BatchDiff(Range&& In);

  void BatchDiff_(Slice In);

  template <typename Range>
  void Flatten(Range&& Out);

  template <typename Range>
  auto KNN(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq);

  auto RangeCount(Box const& query_box);

  template <typename Range>
  auto RangeQuery(Box const& query_box, Range&& Out);

  constexpr static char const* GetTreeName() { return "PTree"; }

  SplitRule split_rule_;

  CpamAugMap cpam_aug_map_;
};

}  // namespace pspt

#include "p_tree_impl/p_batch_delete.hpp"
#include "p_tree_impl/p_batch_diff.hpp"
#include "p_tree_impl/p_batch_insert.hpp"
#include "p_tree_impl/p_build_tree.hpp"
#include "p_tree_impl/p_inter_node.hpp"
#include "p_tree_impl/p_override.hpp"

#endif  // PSPT_P_TREE_H
