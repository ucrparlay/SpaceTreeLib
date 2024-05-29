#include <type_traits>
#include "../oct_tree.h"

namespace cpdd {

template<typename Point>
void octTree<Point>::DeleteTree() {
    this->template DeleteTreeWrapper<leaf<ZValueType, ZValueSlice,std::true_type>, interior>();
}

}  // namespace cpdd
