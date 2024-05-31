#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point>
class KdTree : public BaseTree<Point> {
   public:
    using BaseTree = BaseTree<Point>;

    using BucketType = BaseTree::BucketType;
    using BallsType = BaseTree::BallsType;
    using DimsType = BaseTree::DimsType;
    using Coord = typename Point::Coord;
    using Coords = typename Point::Coords;
    using AugType = bool;
    using Num = Num_Comparator<Coord>;
    using slice = BaseTree::slice;
    using Points = BaseTree::Points;
    using PointsIter = BaseTree::PointsIter;
    using Splitter = BaseTree::Splitter;
    using SplitterSeq = BaseTree::SplitterSeq;
    using Box = BaseTree::Box;
    using BoxSeq = BaseTree::BoxSeq;
    using Circle = BaseTree::Circle;

    struct interior;
    // NOTE: wrapper
    void build(slice In, const DimsType DIM) override;

    void deleteTree() override;
};

}  // namespace cpdd

#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"
