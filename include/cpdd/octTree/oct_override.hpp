#include "../oct_tree.h"

namespace cpdd {

template<typename point>
void
octTree<point>::delete_tree() {
  this->template delete_tree_wrapper<interior>();
}

}  // namespace cpdd
