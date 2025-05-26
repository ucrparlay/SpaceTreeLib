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
  // using AugType = std::optional<bool>;
  using SplitRuleType = SplitRule;

  // NOTE: for the CPAM
  using CurveCode = typename Point::AT::CurveCode;
  using IdType = typename Point::AT::IdType;

  using CpamKey = typename Point::AT;  // morton_id, id
  using CpamVal = Coords;
  // using CpamAug = std::pair<Box, size_t>;
  // PERF: aug value can remove the size, because the tree has this info
  using CpamAug = Box;

  struct CpamEntry {
    using key_t = CpamKey;
    using val_t = CpamVal;
    using aug_t = CpamAug;
    using entry_t = Point;

    using filling_curve_t = SplitRule;
    using sort_output_value_t = std::pair<typename Point::AT, entry_t*>;
    // using entry_t_ref_v = Point*;
    using entry_t_ref_wrapper_v = std::reference_wrapper<Point>;

    static inline key_t const& get_key(entry_t const& e) { return e.GetAug(); }
    static inline val_t const& get_val(entry_t const& e) {
      return e.GetCoords();
    }
    static inline void set_val(entry_t& e, val_t const& v) {
      e.GetCoords() = v;
    }
    static inline entry_t to_entry(key_t const& k, val_t const& v) {
      return entry_t(v, k);
      // return std::make_tuple(k, v);
    };
    static inline aug_t from_entry(entry_t const& e) {
      return from_entry(get_key(e), get_val(e));
    }

    // how to compare key
    static inline bool comp(key_t const& a, key_t const& b) { return a < b; }

    // get an empty aug
    static aug_t get_empty() { return BT::GetEmptyBox(); }

    // WARN: this invoke implicity conversion from Coords to BasicPoint
    static aug_t from_entry(key_t const& k, val_t const& v) {
      return Box(v, v);
    }

    // combine two aug val to get a new aug
    static aug_t combine(aug_t const& a, aug_t const& b) {
      return BT::GetBox(a, b);
    }
  };

  using CpamAugMap = cpam::aug_map<CpamEntry, 16>;
  using CpamMap = CpamAugMap::Map;

  using CpamInnerEntryType =
      std::tuple<typename CpamEntry::key_t,
                 std::reference_wrapper<typename CpamEntry::val_t>>;

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

  // TODO: need to hide the basetree assests, e.g., root_

  SplitRule space_filling_curve_;

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
