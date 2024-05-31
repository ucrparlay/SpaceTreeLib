#pragma once

#include <parlay/slice.h>
#include <functional>
#include <utility>
#include <algorithm>
#include "../kd_tree.h"
#include "parlay/primitives.h"
#include "typeinfo"

namespace cpdd {
template<typename Point>
void KdTree<Point>::build(slice A, const DimsType DIM) {
    // build_z_value(A, DIM);
    // build_z_value_pointer(A, DIM);
    // build_point_z_value(A, DIM);
    // build_point(A, DIM);
    return;
}

}  // namespace cpdd
