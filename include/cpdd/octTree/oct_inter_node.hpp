#pragma once

#include "../oct_tree.h"

namespace cpdd {

template<typename point>
template<typename aug_type>
struct octTree<point>::interior : node {
    node* left;
    node* right;
    z_bit_type bit;
    aug_type aug;
    interior(node* _left, node* _right, z_bit_type _bit) :
        node{false, _left->size + _right->size}, left(_left), right(_right), bit(_bit), aug(false) {}
};

template<typename point>
template<typename aug_type>
typename octTree<point>::template interior<aug_type>* octTree<point>::alloc_oct_interior_node(node* L, node* R,
                                                                                              z_bit_type bit) {
    interior<aug_type>* o = parlay::type_allocator<interior<aug_type>>::alloc();
    new (o) interior(L, R, bit);
    return o;
}

}  // namespace cpdd
