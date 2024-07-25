#include <type_traits>
#include "../kd_tree.h"

namespace cpdd {

template<typename Point>
void KdTree<Point>::DeleteTree() {
    this->template DeleteTreeWrapper<Leaf<Point, Slice, std::true_type>, Interior<Point, bool>>();
}

}  // namespace cpdd
