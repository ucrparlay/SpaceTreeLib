#pragma once

#include <parlay/parallel.h>
#include "../base_tree.h"

namespace cpdd {

template<typename point>
template<typename interior>
bool baseTree<point>::checkBox(node* T, const box& bx) {
  assert(T != nullptr);
  assert(legal_box(bx));
  points wx = points::uninitialized(T->size);
  flatten(T, parlay::make_slice(wx));
  auto b = get_box(parlay::make_slice(wx));
  // LOG << b.first << b.second << ENDL;
  return within_box(b, bx);
}

template<typename point>
template<typename interior>
size_t baseTree<point>::checkSize(node* T) {
  if (T->is_leaf) {
    return T->size;
  }
  interior* TI = static_cast<interior*>(T);
  size_t l = checkSize<interior>(TI->left);
  size_t r = checkSize<interior>(TI->right);
  assert(l + r == T->size);
  return T->size;
}

template<typename point>
template<typename interior>
void baseTree<point>::checkTreeSameSequential(node* T, int dim,
                                              const int& DIM) {
  if (T->is_leaf) {
    // assert( pick_rebuild_dim( T, DIM ) == dim );
    return;
  }
  interior* TI = static_cast<interior*>(T);
  if (TI->split.second != dim) {
    LOG << int(TI->split.second) << " " << int(dim) << TI->size << ENDL;
  }
  assert(TI->split.second == dim);
  dim = (dim + 1) % DIM;
  parlay::par_do_if(
      T->size > 1000,
      [&]() { checkTreeSameSequential<interior>(TI->left, dim, DIM); },
      [&]() { checkTreeSameSequential<interior>(TI->right, dim, DIM); });
  return;
}

template<typename point>
template<typename interior>
void baseTree<point>::validate(const dim_type DIM) {
  if (checkBox<interior>(this->root, this->bbox) && legal_box(this->bbox)) {
    std::cout << "Correct bounding box" << std::endl << std::flush;
  } else {
    std::cout << "wrong bounding box" << std::endl << std::flush;
    abort();
  }

  if (this->_split_rule == ROTATE_DIM) {
    checkTreeSameSequential<interior>(this->root, 0, DIM);
    std::cout << "Correct rotate dimension" << std::endl << std::flush;
  }

  if (checkSize<interior>(this->root) == this->root->size) {
    std::cout << "Correct size" << std::endl << std::flush;
  } else {
    std::cout << "wrong tree size" << std::endl << std::flush;
    abort();
  }
  return;
}

template<typename point>
template<typename interior>
size_t baseTree<point>::getTreeHeight() {
  size_t deep = 0;
  return getMaxTreeDepth<interior>(this->root, deep);
}

template<typename point>
template<typename interior>
size_t baseTree<point>::getMaxTreeDepth(node* T, size_t deep) {
  if (T->is_leaf) {
    return deep;
  }
  interior* TI = static_cast<interior*>(T);
  int l = getMaxTreeDepth<interior>(TI->left, deep + 1);
  int r = getMaxTreeDepth<interior>(TI->right, deep + 1);
  return std::max(l, r);
}

template<typename point>
template<typename interior>
double baseTree<point>::getAveTreeHeight() {
  parlay::sequence<size_t> heights(this->root->size);
  size_t idx = 0;
  countTreeHeights<interior>(this->root, 0, idx, heights);
  // auto kv = parlay::histogram_by_key( heights.cut( 0, idx ) );
  // std::sort( kv.begin(), kv.end(),
  //            [&]( auto a, auto b ) { return a.first < b.first; } );
  // for ( auto i : kv )
  //     LOG << i.first << " " << i.second << ENDL;
  return double(1.0 * parlay::reduce(heights.cut(0, idx)) / idx);
}

template<typename point>
template<typename interior>
size_t baseTree<point>::countTreeNodesNum(node* T) {
  if (T->is_leaf) {
    return 1;
  }

  interior* TI = static_cast<interior*>(T);
  size_t l, r;
  parlay::par_do([&]() { l = countTreeNodesNum<interior>(TI->left); },
                 [&]() { r = countTreeNodesNum<interior>(TI->right); });
  return l + r + 1;
}

template<typename point>
template<typename interior>
void baseTree<point>::countTreeHeights(node* T, size_t deep, size_t& idx,
                                       parlay::sequence<size_t>& heights) {
  if (T->is_leaf) {
    heights[idx++] = deep;
    return;
  }
  interior* TI = static_cast<interior*>(T);
  countTreeHeights<interior>(TI->left, deep + 1, idx, heights);
  countTreeHeights<interior>(TI->right, deep + 1, idx, heights);
  return;
}

}  // namespace cpdd
