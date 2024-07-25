#pragma once

#include <parlay/range.h>
#include <parlay/type_traits.h>
#include <parlay/slice.h>
#include "../kd_tree.h"

namespace cpdd {
template<typename Point>
template<typename R>
void KdTree<Point>::Build(R&& In, int DIM) {
    static_assert(parlay::is_random_access_range_v<R>);
    static_assert(parlay::is_less_than_comparable_v<parlay::range_reference_type_t<R>>);
    static_assert(std::is_constructible_v<parlay::range_value_type_t<R>, parlay::range_reference_type_t<R>>);
    Slice A = parlay::make_slice(In);
    Build_(A, DIM);
}

template<typename Point>
void KdTree<Point>::Build_(Slice A, const DimsType DIM) {
    Points B = Points::uninitialized(A.size());
    this->tree_box_ = BT::GetBox(A);
    // this->root = build_recursive(A, B.cut(0, A.size()), 0, DIM, this->bbox);
    assert(this->root_ != nullptr);
    return;
}

}  // namespace cpdd
