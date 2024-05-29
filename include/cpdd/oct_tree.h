#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"
#include "libmorton/morton.h"

namespace cpdd {

template<typename Point>
class octTree : public BaseTree<Point> {
   public:
    using BaseTree = BaseTree<Point>;

    using BucketType = BaseTree::BucketType;
    using BallsType = BaseTree::BallsType;
    using DimsType = BaseTree::DimsType;
    using ZBitType = uint_fast8_t;
    using ZValueType = uint_fast64_t;
    using ZValueSlice = parlay::slice<ZValueType*, ZValueType*>;
    using ZValuePointerPair = std::pair<ZValueType, Point*>;
    using ZValuePointerSlice = parlay::slice<ZValuePointerPair*, ZValuePointerPair*>;
    using ZValuePointPair = std::pair<ZValueType, Point>;
    using ZValuePointSlice = parlay::slice<ZValuePointPair*, ZValuePointPair*>;

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

    static constexpr ZBitType KEY_BITS = 64;

    struct interior;

    static interior* alloc_oct_interior_node(node* L, node* R, ZBitType bit);

    static inline ZValueType interleave_bits(Point* p, const DimsType DIM);

    static inline ZValueType get_z_value(const Point& p);

    void build_z_value(slice In, const DimsType DIM);
    node* build_recursive_with_z_value(ZValueSlice In, ZBitType bit, const DimsType DIM);

    void build_z_value_pointer(slice In, const DimsType DIM);
    node* build_recursive_with_z_value_pointer(ZValuePointerSlice In, ZBitType bit, const DimsType DIM);

    void build_point_z_value(slice In, const DimsType DIM);
    node* build_recursive_point_z_value(ZValuePointSlice In, ZBitType bit, const DimsType DIM);

    void build_point(slice In, const DimsType DIM);
    node* build_recursive_point(slice In, ZBitType bit, const DimsType DIM);

    // NOTE: wrapper
    void Build(slice In, const DimsType DIM) override;
    node* serial_build_recursive(slice In, ZBitType bit, const DimsType DIM);
    node* build_recursive(slice In, ZBitType bit, const DimsType DIM);

    void DeleteTree() override;

    uint64_t binary_search_time = 0;
    uint64_t leaf_alloc_time = 0;
    uint64_t leaf_num = 0;
    uint64_t time_base = 1000000;
};

}  // namespace cpdd

#include "oct_tree_impl/oct_build_tree.hpp"
#include "oct_tree_impl/oct_inter_node.hpp"
#include "oct_tree_impl/oct_override.hpp"
