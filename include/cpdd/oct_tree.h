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
    using z_bit_type = uint_fast8_t;
    using z_value_type = uint_fast64_t;
    using z_value_slice = parlay::slice<z_value_type*, z_value_type*>;
    using z_value_pointer_pair = std::pair<z_value_type, Point*>;
    using z_value_pointer_slice = parlay::slice<z_value_pointer_pair*, z_value_pointer_pair*>;
    using z_value_point_pair = std::pair<z_value_type, Point>;
    using z_value_point_slice = parlay::slice<z_value_point_pair*, z_value_point_pair*>;

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

    static constexpr z_bit_type KEY_BITS = 64;

    struct interior;

    static interior* alloc_oct_interior_node(node* L, node* R, z_bit_type bit);

    static inline z_value_type interleave_bits(Point* p, const DimsType DIM);

    static inline z_value_type get_z_value(const Point& p);

    void build_z_value(slice In, const DimsType DIM);
    node* build_recursive_with_z_value(z_value_slice In, z_bit_type bit, const DimsType DIM);

    void build_z_value_pointer(slice In, const DimsType DIM);
    node* build_recursive_with_z_value_pointer(z_value_pointer_slice In, z_bit_type bit, const DimsType DIM);

    void build_point_z_value(slice In, const DimsType DIM);
    node* build_recursive_point_z_value(z_value_point_slice In, z_bit_type bit, const DimsType DIM);

    void build_point(slice In, const DimsType DIM);
    node* build_recursive_point(slice In, z_bit_type bit, const DimsType DIM);

    // NOTE: wrapper
    void Build(slice In, const DimsType DIM) override;
    node* serial_build_recursive(slice In, z_bit_type bit, const DimsType DIM);
    node* build_recursive(slice In, z_bit_type bit, const DimsType DIM);

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
