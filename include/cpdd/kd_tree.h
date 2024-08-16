#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint_fast8_t kBDO = 6>
class KdTree : private BaseTree<Point, kBDO> {
 public:
    using BT = BaseTree<Point, kBDO>;

    using BucketType = BT::BucketType;
    using BallsType = BT::BallsType;
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
    using NodeTag = BT::NodeTag;
    using TagNodes = BT::TagNodes;
    using NodeBox = BT::NodeBox;
    using Splitter = HyperPlane;
    using SplitterSeq = HyperPlaneSeq;

    using AugType = bool;
    struct KdInteriorNode;

    using SplitRuleType = SplitRule;
    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = KdInteriorNode;

    // NOTE: expose basetree interface
    using BT::GetAveTreeHeight;
    using BT::GetBox;
    using BT::GetMaxTreeDepth;
    using BT::GetRoot;
    using BT::GetRootBox;
    using BT::SetRoot;
    using BT::Validate;

    // NOTE: functions
    template<typename Range>
    void Build(Range&& In);

    void DeleteTree() override;

    void BatchInsert(Slice In, const DimsType BT::kDim);

    Node* RebuildWithInsert(Node* T, Slice In, const DimsType d,
                            const DimsType BT::kDim);

    static inline void UpdateInterior(Node* T, Node* L, Node* R);

    Node* BatchInsertRecursive(Node* T, Slice In, Slice Out, DimsType d,
                               const DimsType BT::kDim);

    static Node* UpdateInnerTreeByTag(BucketType idx, const NodeTag& tags,
                                      parlay::sequence<Node*>& treeNodes,
                                      BucketType& p, const TagNodes& rev_tag);

    // NOTE: batch delete
    // NOTE: in default, all Points to be deleted are assumed in the tree
    void BatchDelete(Slice In, const DimsType BT::kDim);

    // NOTE: explicitly specify all Points to be deleted are in the tree
    void BatchDelete(Slice In, const DimsType BT::kDim, FullCoveredTag);

    // NOTE: for the case that some Points to be deleted are not in the tree
    void BatchDelete(Slice In, const DimsType BT::kDim, PartialCoverTag);

    //  PERF: try pass a reference to bx
    NodeBox BatchDeleteRecursive(Node* T, const Box& bx, Slice In, Slice Out,
                                  DimsType d, const DimsType BT::kDim, bool hasTomb,
                                  FullCoveredTag);

    // TODO: add bounding Box for batch delete recursive as well
    // WARN: fix the possible in partial deletion as well
    NodeBox BatchDeleteRecursive(Node* T, const Box& bx, Slice In, Slice Out,
                                  DimsType d, const DimsType BT::kDim,
                                  PartialCoverTag);

    NodeBox DeleteInnerTree(BucketType idx, const NodeTag& tags,
                              parlay::sequence<NodeBox>& treeNodes,
                              BucketType& p, const TagNodes& rev_tag,
                              const DimsType d, const DimsType BT::kDim);

    template<typename Range>
    void Flatten(Range&& Out);

    template<typename Range>
    auto KNN(Node* T, const Point& q, kBoundedQueue<Point, Range>& bq);

    size_t RangeCount(const Box& query_box);

    size_t RangeCount(const Circle& cl);

    template<typename Range>
    size_t RangeQuery(const Box& query_box, Range&& Out);

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
