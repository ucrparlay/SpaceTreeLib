#include <type_traits>
#include "../oct_tree.h"

namespace cpdd {

template<typename point>
void octTree<point>::delete_tree() {
  this->template delete_tree_wrapper<
      leaf<z_value_type, z_value_slice, std::true_type>, interior>();
}

}  // namespace cpdd
