#pragma once

#include "../oct_tree.h"

namespace cpdd {

template<typename Point>
struct octTree<Point>::interior : Node {
    Node* left;
    Node* right;
    ZBitType bit;
    AugType aug;
    interior(Node* _left, Node* _right, ZBitType _bit) :
        Node{false, _left->size + _right->size}, left(_left), right(_right), bit(_bit), aug(false) {}
};

template<typename Point>
typename octTree<Point>::interior* octTree<Point>::alloc_oct_interior_node(Node* L, Node* R, ZBitType bit) {
    interior* o = parlay::type_allocator<interior>::alloc();
    new (o) interior(L, R, bit);
    return o;
}

}  // namespace cpdd
