#pragma once

#include <functional>
#include <utility>
#include <optional>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint_fast8_t kBDO = 6>
class KdTree : private BaseTree<Point, KdTree<Point, SplitRule, kBDO>, kBDO> {
 public:
    using BT = BaseTree<Point, KdTree<Point, SplitRule, kBDO>, kBDO>;

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
    using Tag2Node = typename BT::Tag2Node;
    using NodeBox = typename BT::NodeBox;
    using NodeBoxSeq = typename BT::NodeBoxSeq;
    using Splitter = HyperPlane;
    using SplitterSeq = HyperPlaneSeq;

    using AugType = std::optional<bool>;
    // using AugType = bool;
    struct KdInteriorNode;

    using SplitRuleType = SplitRule;
    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = KdInteriorNode;
    using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
    using BoxCut = typename BT::BoxCut;

    // NOTE: expose basetree interface
    using BT::GetAveTreeHeight;
    using BT::GetBox;
    using BT::GetMaxTreeDepth;
    using BT::GetRoot;
    using BT::GetRootBox;
    using BT::SetRoot;
    using BT::Validate;

    template<typename Leaf, typename Interior, bool granularity,
             typename... Args>
    friend Node* BT::RebuildSingleTree(Node* T, Args&&... args);

    template<typename Leaf, typename Interior, typename... Args>
    friend Node* BT::RebuildWithInsert(Node* T, Slice In, Args&&... args);

    void KdTreeTag();

    // NOTE: functions
    template<typename Range>
    void Build(Range&& In);

    void DeleteTree() override;

    // TODO: move it to inner tree
    static NodeBox UpdateInnerTreePointerBox(BucketType idx,
                                             const NodeTagSeq& tags,
                                             NodeBoxSeq& tree_nodes,
                                             BucketType& p);

    // TODO: move it to inner tree
    static Node* UpdateInnerTreePointer(BucketType idx, const NodeTagSeq& tags,
                                        parlay::sequence<Node*>& tree_nodes,
                                        BucketType& p);

    NodeBox RebuildTreeRecursive(Node* T, DimsType d,
                                 const bool granularity = true);

    void BatchInsert(Slice In);

    Node* BatchInsertRecursive(Node* T, Slice In, Slice Out, DimsType d);

    // NOTE: batch delete
    // NOTE: in default, all Points to be deleted are assumed in the tree
    template<typename Range>
    void BatchDelete(Range&& In);

    void BatchDelete_(Slice In);

    NodeBox BatchDeleteRecursive(Node* T, const Box& bx, Slice In, Slice Out,
                                 DimsType d, bool has_tomb);

    // NOTE: batch diff
    // NOTE: for the case that some Points to be deleted are not in the tree
    template<typename Range>
    void BatchDiff(Range&& In);

    void BatchDiff_(Slice In);

    // TODO: add bounding Box for batch delete recursive as well
    // WARN: fix the possible in partial deletion as well
    NodeBox BatchDiffRecursive(Node* T, const Box& bx, Slice In, Slice Out,
                               DimsType d);

    NodeBox DeleteInnerTree(BucketType idx, const NodeTagSeq& tags,
                            NodeBoxSeq& tree_nodes, BucketType& p,
                            const DimsType d);

    template<typename Range>
    void Flatten(Range&& Out);

    template<typename Range>
    auto KNN(Node* T, const Point& q, kBoundedQueue<Point, Range>& bq);

    auto RangeCount(const Box& query_box);

    auto RangeCount(const Circle& cl);

    template<typename Range>
    auto RangeQuery(const Box& query_box, Range&& Out);

 private:
    void DivideRotate(Slice In, SplitterSeq& pivots, DimsType dim,
                      BucketType idx, BoxSeq& boxs, const Box& bx);

    void PickPivots(Slice In, const size_t& n, SplitterSeq& pivots,
                    const DimsType dim, BoxSeq& boxs, const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, DimsType dim, const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                               const Box& bx);

    void Build_(Slice In);

    SplitRule split_rule_;
};

}  // namespace cpdd

#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"
#include "kd_tree_impl/kd_batch_insert.hpp"
#include "kd_tree_impl/kd_batch_delete.hpp"
#include "kd_tree_impl/kd_batch_diff.hpp"
