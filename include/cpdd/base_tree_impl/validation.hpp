#pragma once

#include <parlay/parallel.h>
#include "../base_tree.h"

namespace cpdd {

template<typename Point>
template<typename Interior>
bool BaseTree<Point>::CheckBox(Node* T, const Box& bx) {
    assert(T != nullptr);
    assert(LegalBox(bx));
    Points wx = Points::uninitialized(T->size);
    Flatten(T, parlay::make_slice(wx));
    auto b = GetBox(parlay::make_slice(wx));
    // LOG << b.first << b.second << ENDL;
    return WithinBox(b, bx);
}

template<typename Point>
template<typename Interior>
size_t BaseTree<Point>::CheckSize(Node* T) {
    if (T->is_leaf) {
        return T->size;
    }
    Interior* TI = static_cast<Interior*>(T);
    size_t l = CheckSize<Interior>(TI->left);
    size_t r = CheckSize<Interior>(TI->right);
    assert(l + r == T->size);
    return T->size;
}

template<typename Point>
template<typename Interior>
void BaseTree<Point>::CheckTreeSameSequential(Node* T, int dim,
                                              const int& DIM) {
    if (T->is_leaf) {
        // assert( PickRebuildDim( T, DIM ) == dim );
        return;
    }
    Interior* TI = static_cast<Interior*>(T);
    if (TI->split.second != dim) {
        LOG << int(TI->split.second) << " " << int(dim) << TI->size << ENDL;
    }
    assert(TI->split.second == dim);
    dim = (dim + 1) % DIM;
    parlay::par_do_if(
        T->size > 1000,
        [&]() { CheckTreeSameSequential<Interior>(TI->left, dim, DIM); },
        [&]() { CheckTreeSameSequential<Interior>(TI->right, dim, DIM); });
    return;
}

template<typename Point>
template<typename Interior>
void BaseTree<Point>::Validate(const DimsType DIM) {
    if (CheckBox<Interior>(this->root_, this->tree_box_) &&
        LegalBox(this->tree_box_)) {
        std::cout << "Correct bounding Box" << std::endl << std::flush;
    } else {
        std::cout << "wrong bounding Box" << std::endl << std::flush;
        abort();
    }

    if (this->split_rule_ == kRotateDim) {
        CheckTreeSameSequential<Interior>(this->root_, 0, DIM);
        std::cout << "Correct rotate dimension" << std::endl << std::flush;
    }

    if (CheckSize<Interior>(this->root_) == this->root_->size) {
        std::cout << "Correct size" << std::endl << std::flush;
    } else {
        std::cout << "wrong tree size" << std::endl << std::flush;
        abort();
    }
    return;
}

template<typename Point>
template<typename Interior>
size_t BaseTree<Point>::GetTreeHeight() {
    size_t deep = 0;
    return GetMaxTreeDepth<Interior>(this->root_, deep);
}

template<typename Point>
template<typename Interior>
size_t BaseTree<Point>::GetMaxTreeDepth(Node* T, size_t deep) {
    if (T->is_leaf) {
        return deep;
    }
    Interior* TI = static_cast<Interior*>(T);
    int l = GetMaxTreeDepth<Interior>(TI->left, deep + 1);
    int r = GetMaxTreeDepth<Interior>(TI->right, deep + 1);
    return std::max(l, r);
}

template<typename Point>
template<typename Interior>
double BaseTree<Point>::GetAveTreeHeight() {
    parlay::sequence<size_t> heights(this->root_->size);
    size_t idx = 0;
    CountTreeHeights<Interior>(this->root_, 0, idx, heights);
    // auto kv = parlay::histogram_by_key( heights.cut( 0, idx ) );
    // std::sort( kv.begin(), kv.end(),
    //            [&]( auto a, auto b ) { return a.first < b.first; } );
    // for ( auto i : kv )
    //     LOG << i.first << " " << i.second << ENDL;
    return double(1.0 * parlay::reduce(heights.cut(0, idx)) / idx);
}

template<typename Point>
template<typename Interior>
size_t BaseTree<Point>::CountTreeNodesNum(Node* T) {
    if (T->is_leaf) {
        return 1;
    }

    Interior* TI = static_cast<Interior*>(T);
    size_t l, r;
    parlay::par_do([&]() { l = CountTreeNodesNum<Interior>(TI->left); },
                   [&]() { r = CountTreeNodesNum<Interior>(TI->right); });
    return l + r + 1;
}

template<typename Point>
template<typename Interior>
void BaseTree<Point>::CountTreeHeights(Node* T, size_t deep, size_t& idx,
                                       parlay::sequence<size_t>& heights) {
    if (T->is_leaf) {
        heights[idx++] = deep;
        return;
    }
    Interior* TI = static_cast<Interior*>(T);
    CountTreeHeights<Interior>(TI->left, deep + 1, idx, heights);
    CountTreeHeights<Interior>(TI->right, deep + 1, idx, heights);
    return;
}

}  // namespace cpdd
