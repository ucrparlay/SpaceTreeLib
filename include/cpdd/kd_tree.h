#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule>
class KdTree : public BaseTree<Point> {
 public:
    using BT = BaseTree<Point>;

    using BucketType = BT::BucketType;
    using BallsType = BT::BallsType;
    using DimsType = BT::DimsType;
    using Coord = typename Point::Coord;
    using Coords = typename Point::Coords;
    using AugType = bool;
    using Num = Num_Comparator<Coord>;
    using Slice = BT::Slice;
    using Points = BT::Points;
    using PointsIter = BT::PointsIter;
    using Splitter = BT::Splitter;
    using SplitterSeq = BT::SplitterSeq;
    using Box = BT::Box;
    using BoxSeq = BT::BoxSeq;
    using Circle = BT::Circle;

    struct interior;

    template<typename Range>
    void Build(Range&& In, int DIM);

    void Build_(Slice In, const DimsType DIM);

    void DeleteTree() override;

 private:
    void DivideRotate(Slice In, SplitterSeq& pivots, DimsType dim,
                      BucketType idx, BucketType deep, BucketType& bucket,
                      const DimsType DIM, BoxSeq& boxs, const Box& bx);

    void PickPivots(Slice In, const size_t& n, SplitterSeq& pivots,
                    const DimsType dim, const DimsType DIM, BoxSeq& boxs,
                    const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, DimsType dim, const DimsType DIM,
                         const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                               const DimsType DIM, const Box& bx);

    SplitRule split_rule_;
};

}  // namespace cpdd

#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"
