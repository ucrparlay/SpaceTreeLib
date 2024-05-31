#include <type_traits>
#include "../oct_tree.h"

namespace cpdd {

template<typename Point>
void octTree<Point>::deleteTree() {
    this->template deleteTreeWrapper<Leaf<ZValueType, ZValueSlice, std::true_type>, interior>();
}

}  // namespace cpdd
