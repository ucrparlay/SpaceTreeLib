#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point>
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

    template<typename R>
    void Build(R&& In, int DIM);

    void Build_(Slice In, const DimsType DIM) override;

    void divide_rotate(Slice In, SplitterSeq& pivots, DimsType dim, BucketType idx, BucketType deep, BucketType& bucket,
                       const DimsType DIM, BoxSeq& boxs, const Box& bx);
    void pick_pivots(Slice In, const size_t& n, SplitterSeq& pivots, const DimsType dim, const DimsType DIM,
                     BoxSeq& boxs, const Box& bx);
    static inline BucketType find_bucket(const Point& p, const SplitterSeq& pivots);
    static void partition(Slice A, Slice B, const size_t n, const SplitterSeq& pivots,
                          parlay::sequence<BallsType>& sums);
    static Node* build_inner_tree(BucketType idx, SplitterSeq& pivots, parlay::sequence<Node*>& treeNodes);
    void build(Slice In, const DimsType DIM);
    PointsIter serial_partition(Slice In, DimsType d);
    Node* serial_build_recursive(Slice In, Slice Out, DimsType dim, const DimsType DIM, const Box& bx);
    Node* build_recursive(Slice In, Slice Out, DimsType dim, const DimsType DIM, const Box& bx);

    void DeleteTree() override;
};

}  // namespace cpdd

#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"
