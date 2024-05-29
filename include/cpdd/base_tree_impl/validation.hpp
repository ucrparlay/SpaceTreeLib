#pragma once

#include <parlay/parallel.h>
#include "../base_tree.h"

namespace cpdd {

template<typename Point>
template<typename interior>
bool BaseTree<Point>::CheckBox(node* T, const Box& bx) {
    assert(T != nullptr);
    assert(LegalBox(bx));
    Points wx = Points::uninitialized(T->size);
    Flatten(T, parlay::make_slice(wx));
    auto b = GetBox(parlay::make_slice(wx));
    // LOG << b.first << b.second << ENDL;
    return WithinBox(b, bx);
}

template<typename Point>
template<typename interior>
size_t BaseTree<Point>::CheckSize(node* T) {
    if (T->is_leaf) {
        return T->size;
    }
    interior* TI = static_cast<interior*>(T);
    size_t l = CheckSize<interior>(TI->left);
    size_t r = CheckSize<interior>(TI->right);
    assert(l + r == T->size);
    return T->size;
}

template<typename Point>
template<typename interior>
void BaseTree<Point>::CheckTreeSameSequential(node* T, int dim, const int& DIM) {
    if (T->is_leaf) {
        // assert( PickRebuildDim( T, DIM ) == dim );
        return;
    }
    interior* TI = static_cast<interior*>(T);
    if (TI->split.second != dim) {
        LOG << int(TI->split.second) << " " << int(dim) << TI->size << ENDL;
    }
    assert(TI->split.second == dim);
    dim = (dim + 1) % DIM;
    parlay::par_do_if(
        T->size > 1000, [&]() { CheckTreeSameSequential<interior>(TI->left, dim, DIM); },
        [&]() { CheckTreeSameSequential<interior>(TI->right, dim, DIM); });
    return;
}

template<typename Point>
template<typename interior>
void BaseTree<Point>::Validate(const DimsType DIM) {
    if (CheckBox<interior>(this->root_, this->tree_box_) && LegalBox(this->tree_box_)) {
        std::cout << "Correct bounding Box" << std::endl << std::flush;
    } else {
        std::cout << "wrong bounding Box" << std::endl << std::flush;
        abort();
    }

    if (this->split_rule_ == kRotateDim) {
        CheckTreeSameSequential<interior>(this->root_, 0, DIM);
        std::cout << "Correct rotate dimension" << std::endl << std::flush;
    }

    if (CheckSize<interior>(this->root_) == this->root_->size) {
        std::cout << "Correct size" << std::endl << std::flush;
    } else {
        std::cout << "wrong tree size" << std::endl << std::flush;
        abort();
    }
    return;
}

template<typename Point>
template<typename interior>
size_t BaseTree<Point>::GetTreeHeight() {
    size_t deep = 0;
    return GetMaxTreeDepth<interior>(this->root_, deep);
}

template<typename Point>
template<typename interior>
size_t BaseTree<Point>::GetMaxTreeDepth(node* T, size_t deep) {
    if (T->is_leaf) {
        return deep;
    }
    interior* TI = static_cast<interior*>(T);
    int l = GetMaxTreeDepth<interior>(TI->left, deep + 1);
    int r = GetMaxTreeDepth<interior>(TI->right, deep + 1);
    return std::max(l, r);
}

template<typename Point>
template<typename interior>
double BaseTree<Point>::GetAveTreeHeight() {
    parlay::sequence<size_t> heights(this->root_->size);
    size_t idx = 0;
    CountTreeHeights<interior>(this->root_, 0, idx, heights);
    // auto kv = parlay::histogram_by_key( heights.cut( 0, idx ) );
    // std::sort( kv.begin(), kv.end(),
    //            [&]( auto a, auto b ) { return a.first < b.first; } );
    // for ( auto i : kv )
    //     LOG << i.first << " " << i.second << ENDL;
    return double(1.0 * parlay::reduce(heights.cut(0, idx)) / idx);
}

template<typename Point>
template<typename interior>
size_t BaseTree<Point>::CountTreeNodesNum(node* T) {
    if (T->is_leaf) {
        return 1;
    }

    interior* TI = static_cast<interior*>(T);
    size_t l, r;
    parlay::par_do([&]() { l = CountTreeNodesNum<interior>(TI->left); },
                   [&]() { r = CountTreeNodesNum<interior>(TI->right); });
    return l + r + 1;
}

template<typename Point>
template<typename interior>
void BaseTree<Point>::CountTreeHeights(node* T, size_t deep, size_t& idx, parlay::sequence<size_t>& heights) {
    if (T->is_leaf) {
        heights[idx++] = deep;
        return;
    }
    interior* TI = static_cast<interior*>(T);
    CountTreeHeights<interior>(TI->left, deep + 1, idx, heights);
    CountTreeHeights<interior>(TI->right, deep + 1, idx, heights);
    return;
}

}  // namespace cpdd
