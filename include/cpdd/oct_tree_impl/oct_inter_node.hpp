#pragma once

#include "../oct_tree.h"

namespace cpdd {

template<typename Point>
struct octTree<Point>::interior : node {
    node* left;
    node* right;
    ZBitType bit;
    AugType aug;
    interior(node* _left, node* _right, ZBitType _bit) :
        node{false, _left->size + _right->size}, left(_left), right(_right), bit(_bit), aug(false) {}
};

template<typename Point>
typename octTree<Point>::interior* octTree<Point>::alloc_oct_interior_node(node* L, node* R, ZBitType bit) {
    interior* o = parlay::type_allocator<interior>::alloc();
    new (o) interior(L, R, bit);
    return o;
}

}  // namespace cpdd
