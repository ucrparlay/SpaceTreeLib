#include "../base_tree.h"

namespace cpdd {

template<typename point>
template<typename interior_node_type>
void
baseTree<point>::delete_tree_wrapper() {
  if (this->root == nullptr) {
    return;
  }
  delete_tree_recursive<interior_node_type>(this->root);
  this->root = nullptr;
  return;
}

template<typename point>  //* delete tree in parallel
template<typename interior_node_type>
void
baseTree<point>::delete_tree_recursive(node* T, bool granularity) {
  if (T == nullptr) return;
  if (T->is_leaf) {
    free_node<point, leaf>(T);
  } else {
    interior_node_type* TI = static_cast<interior_node_type*>(T);
    // NOTE: enable granularity control by default, if it is disabled, always
    // delete in parallel
    parlay::par_do_if(
        (granularity && T->size > SERIAL_BUILD_CUTOFF) || !granularity,
        [&] {
          delete_tree_recursive<interior_node_type>(TI->left, granularity);
        },
        [&] {
          delete_tree_recursive<interior_node_type>(TI->right, granularity);
        });
    free_node<point, interior_node_type>(T);
  }
}

}  // namespace cpdd
