#pragma once

#include <parlay/slice.h>
#include <functional>
#include <utility>
#include <algorithm>
#include "../oct_tree.h"
#include "libmorton/morton_BMI.h"
#include "parlay/primitives.h"
#include "typeinfo"

namespace cpdd {
template<typename Point>
void octTree<Point>::build(slice A, const DimsType DIM) {
    // build_z_value(A, DIM);
    // build_z_value_pointer(A, DIM);
    // build_point_z_value(A, DIM);
    // build_point(A, DIM);
    return;
}

template<typename Point>
typename octTree<Point>::ZValueType octTree<Point>::interleave_bits(Point* p, const DimsType DIM) {
    ZBitType loc = 0;
    ZValueType id = 0;
    for (ZBitType i = 0; i < KEY_BITS / DIM; i++) {
        for (DimsType d = 0; d < DIM; d++) {
            id = id | (((p->pnt[d] >> i) & static_cast<ZValueType>(1)) << (loc++));
        }
    }
    return id;
}

template<typename Point>
typename octTree<Point>::ZValueType octTree<Point>::get_z_value(const Point& p) {
    if constexpr (std::tuple_size<typename Point::Coords>::value == 2) {
        return libmorton::m2D_e_BMI<ZValueType, Coord>(p.pnt[0], p.pnt[1]);
    } else {
        return libmorton::m3D_e_BMI<ZValueType, Coord>(p.pnt[0], p.pnt[1], p.pnt[2]);
    }
}

// NOTE: 2: z value and the pointer
template<typename Point>
void octTree<Point>::build_z_value_pointer(slice A, const DimsType DIM) {
    parlay::internal::timer t;
    LOG << "here" << ENDL;
    t.start();
    this->tree_box_ = this->GetBox(A);
    t.next("GetBox");

    auto z_value_arr = parlay::map(A, [&](Point& p) { return std::make_pair(this->get_z_value(p), &p); });
    assert(z_value_arr[0].first == this->interleave_bits(&A[0], DIM));
    t.next("generate z_value_pointer");

    parlay::internal::integer_sort_inplace(parlay::make_slice(z_value_arr), [&](const auto& val) { return val.first; });
    t.next("integer_sort_inplace");

    this->root_ = build_recursive_with_z_value_pointer(parlay::make_slice(z_value_arr), DIM * (KEY_BITS / DIM), DIM);
    t.next("build_recursive_z_value");

    assert(this->root_ != nullptr);
    return;
}

//  NOTE: 2.1: z value and the pointer
template<typename Point>
Node* octTree<Point>::build_recursive_with_z_value_pointer(ZValuePointerSlice In, ZBitType bit, const DimsType DIM) {
    size_t n = In.size();

    if (n <= BaseTree::kLeaveWrap) {
        return allocLeafNode<Point, ZValuePointerSlice, AllocNormalLeafTag, parlay::uninitialized_relocate_tag>(
            In, std::max(In.size(), static_cast<size_t>(BaseTree::kLeaveWrap)));
    }

    if (bit == 0) {
        return allocLeafNode<Point, ZValuePointerSlice, alloc_fat_leaf_tag, parlay::uninitialized_relocate_tag>(
            In, In.size());
    }

    ZValueType val = (static_cast<ZValueType>(1)) << (bit - 1);
    ZValueType mask = (bit == 64) ? ~(static_cast<ZValueType>(0)) : ~(~(static_cast<ZValueType>(0)) << bit);
    auto less = [&](const ZValuePointerPair& x) { return (x.first & mask) < val; };

    size_t pos = parlay::internal::binary_search(In, less);

    if (pos == 0 || pos == n) {
        return build_recursive_with_z_value_pointer(In, bit - 1, DIM);
    }

    Node *L, *R;
    parlay::par_do_if(
        n >= this->kSerialBuildCutoff, [&] { L = build_recursive_with_z_value_pointer(In.cut(0, pos), bit - 1, DIM); },
        [&] { R = build_recursive_with_z_value_pointer(In.cut(pos, n), bit - 1, DIM); });
    return alloc_oct_interior_node(L, R, bit);
}
// NOTE: 1: z value only
template<typename Point>
void octTree<Point>::build_z_value(slice A, const DimsType DIM) {
    parlay::internal::timer t;
    LOG << "here" << ENDL;
    t.start();
    this->tree_box_ = this->GetBox(A);
    t.next("GetBox");

    auto z_value_arr = parlay::map(A, [&](const Point& p) { return this->get_z_value(p); });
    t.next("generate_z_value");

    parlay::internal::integer_sort_inplace(parlay::make_slice(z_value_arr), [&](const ZValueType& val) { return val; });
    t.next("integer_sort_inplace");

    this->root_ = build_recursive_with_z_value(parlay::make_slice(z_value_arr), DIM * (KEY_BITS / DIM), DIM);
    t.next("build_recursive_z_value");

    assert(this->root_ != nullptr);
    return;
}

template<typename Point>
Node* octTree<Point>::build_recursive_with_z_value(ZValueSlice In, ZBitType bit, const DimsType DIM) {
    size_t n = In.size();
    if (bit == 0 || n <= BaseTree::kLeaveWrap) {
        // BUG:: need to check the type of the leaf
        return allocLeafNode<ZValueType, ZValueSlice, AllocNormalLeafTag, parlay::uninitialized_copy_tag>(
            In, std::max(In.size(), static_cast<size_t>(BaseTree::kLeaveWrap)));
    }

    ZValueType val = (static_cast<ZValueType>(1)) << (bit - 1);
    ZValueType mask = (bit == 64) ? ~(static_cast<ZValueType>(0)) : ~(~(static_cast<ZValueType>(0)) << bit);
    auto less = [&](const ZValueType& x) { return (x & mask) < val; };
    size_t pos = parlay::internal::binary_search(In, less);

    if (pos == 0 || pos == n) {
        return build_recursive_with_z_value(In, bit - 1, DIM);
    }

    Node *L, *R;
    parlay::par_do_if(
        n >= this->kSerialBuildCutoff, [&] { L = build_recursive_with_z_value(In.cut(0, pos), bit - 1, DIM); },
        [&] { R = build_recursive_with_z_value(In.cut(pos, n), bit - 1, DIM); });
    return alloc_oct_interior_node(L, R, bit);
}

// NOTE: 3: sort with z-value
template<typename Point>
void octTree<Point>::build_point_z_value(slice A, const DimsType DIM) {
    parlay::internal::timer t;
    LOG << "here" << ENDL;
    t.start();
    this->tree_box_ = this->GetBox(A);
    t.next("GetBox");

    auto z_value_arr = parlay::map(A, [&](Point& p) { return std::make_pair(this->get_z_value(p), p); });
    t.next("generate z_value_pointer");

    parlay::internal::integer_sort_inplace(parlay::make_slice(z_value_arr), [&](const auto& val) { return val.first; });
    // parlay::sort_inplace(z_value_arr);
    t.next("integer_sort_inplace");

    this->root_ = build_recursive_point_z_value(parlay::make_slice(z_value_arr), DIM * (KEY_BITS / DIM), DIM);
    t.next("build_recursive");
    assert(this->root_ != nullptr);
    return;
}

template<typename Point>
Node* octTree<Point>::build_recursive_point_z_value(ZValuePointSlice In, ZBitType bit, const DimsType DIM) {
    size_t n = In.size();
    if (bit == 0 || n <= BaseTree::kLeaveWrap) {
        return allocLeafNode<Point, ZValuePointSlice>(In,
                                                      std::max(In.size(), static_cast<size_t>(BaseTree::kLeaveWrap)));
    }

    ZValueType val = (static_cast<ZValueType>(1)) << (bit - 1);
    ZValueType mask = (bit == 64) ? ~(static_cast<ZValueType>(0)) : ~(~(static_cast<ZValueType>(0)) << bit);
    auto less = [&](const ZValuePointPair& x) { return (x.first & mask) < val; };
    size_t pos = parlay::internal::binary_search(In, less);

    if (pos == 0 || pos == n) {
        return build_recursive_point_z_value(In, bit - 1, DIM);
    }

    Node *L, *R;
    parlay::par_do_if(
        n >= this->kSerialBuildCutoff, [&] { L = build_recursive_point_z_value(In.cut(0, pos), bit - 1, DIM); },
        [&] { R = build_recursive_point_z_value(In.cut(pos, n), bit - 1, DIM); });
    return alloc_oct_interior_node(L, R, bit);
}

// NOTE: 4: sort without explicitly computing the z-value
template<typename Point>
void octTree<Point>::build_point(slice A, const DimsType DIM) {
    parlay::internal::timer t;
    LOG << "here" << ENDL;
    t.start();
    this->tree_box_ = this->GetBox(A);
    t.next("GetBox");

    t.next("generate z_value_pointer");
    // compares the interleaved bits of Points p and q without explicitly
    // interleaving them.  From Timothy Chan.
    auto less = [&](const Point& p, const Point& q) {
        DimsType j, k;
        Coord y, x = 0;
        auto less_msb = [](Coord x, Coord y) { return x < y && x < (x ^ y); };
        for (j = k = 0; k < DIM; ++k)
            if (less_msb(x, y = p.pnt[k] ^ q.pnt[k])) {
                j = k;
                x = y;
            }
        return p.pnt[j] < q.pnt[j];
    };
    parlay::sort_inplace(A, less);
    t.next("integer_sort_inplace");

    this->root_ = build_recursive_point(A, DIM * sizeof(Coord) * 8 - 1, DIM);
    t.next("build_recursive");
    assert(this->root_ != nullptr);
    return;
}

template<typename Point>
Node* octTree<Point>::build_recursive_point(slice In, ZBitType bit, const DimsType DIM) {
    size_t n = In.size();
    if (bit == 0 || n <= BaseTree::kLeaveWrap) {
        return allocLeafNode<Point, slice>(In, std::max(In.size(), static_cast<size_t>(BaseTree::kLeaveWrap)));
    }

    auto bits =
        parlay::delayed::map(In, [&](const Point& p) { return 1 == ((p.pnt[DIM - bit % DIM - 1] >> bit / DIM) & 1); });
    size_t pos = std::lower_bound(bits.begin(), bits.end(), 1) - bits.begin();

    if (pos == 0 || pos == n) {
        return build_recursive_point(In, bit - 1, DIM);
    }

    Node *L, *R;
    parlay::par_do_if(
        n >= this->kSerialBuildCutoff, [&] { L = build_recursive_point(In.cut(0, pos), bit - 1, DIM); },
        [&] { R = build_recursive_point(In.cut(pos, n), bit - 1, DIM); });
    return alloc_oct_interior_node(L, R, bit);
}

}  // namespace cpdd
