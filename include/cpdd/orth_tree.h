#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint_fast8_t kMD = 2,
         uint_fast8_t kBDO = 6>
class OrthTree :
    private BaseTree<Point, OrthTree<Point, SplitRule, kMD, kBDO>, kBDO> {
 public:
    static constexpr size_t kSplitterNum = kMD;
    static constexpr size_t kNodeRegions = 1 << kMD;

    using BT = BaseTree<Point, OrthTree<Point, SplitRule, kMD, kBDO>, kBDO>;

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
    using Tag2Node = BT::Tag2Node;
    using AugType = std::optional<bool>;

    // struct KdInteriorNode;
    struct OrthInteriorNode;

    using SplitRuleType = SplitRule;
    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = OrthInteriorNode;
    using OrthNodeArr = Interior::OrthNodeArr;
    using InnerTree = typename BT::template InnerTree<Leaf, Interior>;

    // NOTE: expose basetree interface
    using BT::Expand2Binary;
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

    void OrthTreeTag();

    // NOTE: functions
    template<typename Range, typename... Args>
    void Build(Range&& In, Args&&... args);

    template<typename Range>
    void BatchInsert(Range&& In);

    void DeleteTree() override;

    // NOTE: batch delete
    // NOTE: in default, all Points to be deleted are assumed in the tree
    template<typename Range>
    void BatchDelete(Range&& In);

    void BatchDelete_(Slice In);

    Node* BatchDeleteRecursive(Node* T, Slice In, Slice Out, const Box& box,
                               bool has_tomb);

    // NOTE: batch diff
    // NOTE: for the case that some Points to be deleted are not in the tree
    template<typename Range>
    void BatchDiff(Range&& In);

    void BatchDiff_(Slice In);

    // TODO: add bounding Box for batch delete recursive as well
    // WARN: fix the possible in partial deletion as well
    NodeBox BatchDiffRecursive(Node* T, const Box& bx, Slice In, Slice Out,
                               DimsType d);

    template<typename Range>
    void Flatten(Range&& Out);

    template<typename Range>
    auto KNN(Node* T, const Point& q, kBoundedQueue<Point, Range>& bq);

    auto RangeCount(const Box& query_box);

    auto RangeCount(const Circle& cl);

    template<typename Range>
    auto RangeQuery(const Box& query_box, Range&& Out);

 private:
    void Build_(Slice In);

    void Build_(Slice In, const Box& box);

    void SerialSplit(Slice In, DimsType dim, DimsType idx, const Box& box,
                     const Splitter& split, parlay::sequence<BallsType>& sums,
                     BoxSeq& box_seq);

    void SerialSplitSkeleton(Node* T, Slice In, DimsType dim, DimsType idx,
                             parlay::sequence<BallsType>& sums);

    void DivideRotate(HyperPlaneSeq& pivots, DimsType dim, BucketType idx,
                      BoxSeq& box_seq, const Box& box);

    void PickPivots(Slice In, const size_t& n, HyperPlaneSeq& pivots,
                    const DimsType dim, BoxSeq& boxs, const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, const Box& bx,
                               bool checked_duplicate);

    static Node* OrthBuildInnerTree(BucketType idx, const HyperPlaneSeq& pivots,
                                    const parlay::sequence<Node*>& tree_nodes);

    void BatchInsert_(Slice In);

    Node* BatchInsertRecursive(Node* T, Slice In, Slice Out);

    static Node* UpdateInnerTreeByTag(BucketType idx, const NodeTagSeq& tags,
                                      parlay::sequence<Node*>& tree_nodes,
                                      BucketType& p);

    SplitRule split_rule_;
    size_t alloc_dummy_num_ = 0;
    size_t alloc_empty_num_ = 0;
    size_t alloc_normal_num_ = 0;
    size_t alloc_interior_num_ = 0;
};

}  // namespace cpdd

#include "orth_tree_impl/orth_build_tree.hpp"
#include "orth_tree_impl/orth_inter_node.hpp"
#include "orth_tree_impl/orth_override.hpp"
#include "orth_tree_impl/orth_batch_insert.hpp"
#include "orth_tree_impl/orth_batch_delete.hpp"
