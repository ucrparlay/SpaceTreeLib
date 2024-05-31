#include <type_traits>
#include "../kd_tree.h"

namespace cpdd {

template<typename Point>
void KdTree<Point>::deleteTree() {
    this->template deleteTreeWrapper<Leaf<Point, slice, std::true_type>, Interior<Point, bool>>();
}

}  // namespace cpdd
