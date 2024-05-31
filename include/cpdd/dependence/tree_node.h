#pragma once
#include <parlay/slice.h>
#include <cstdint>
#include "basic_point.h"
#include "utility.h"

namespace cpdd {

struct Node {
    Node() : is_leaf{false}, size{0} {};
    Node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

    bool is_leaf;
    size_t size;
};

template<typename Point, typename Range, typename PointAssignTag>
struct Leaf : Node {
    using Points = parlay::sequence<Point>;

    // NOTE: default allocator
    Leaf() : Node{true, static_cast<size_t>(0)}, is_dummy(false) {};

    // NOTE: alloc a normal leaf
    Leaf(Range In, const auto alloc_size, AllocNormalLeafTag) :
        Node{true, static_cast<size_t>(In.size())}, is_dummy(false), pts(Points::uninitialized(alloc_size)) {
        assert(In.size() <= alloc_size);
        if constexpr (is_pair<typename Range::value_type>::value) {  // NOTE: if input is pair, assign using the second
                                                                     // value
            static_assert(std::is_pointer<typename Range::value_type::second_type>::value);
            for (int i = 0; i < In.size(); i++) {
                parlay::assign_dispatch(pts[i], *(In[i].second), PointAssignTag());
            }
        } else {  //  NOTE: input is a Point with coordinates
            for (int i = 0; i < In.size(); i++) {
                parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
            }
        }
    }

    // NOTE: alloc a dummy leaf
    Leaf(Range In, AllocDummyLeafTag) :
        Node{true, static_cast<size_t>(In.size())}, is_dummy(true), pts(Points::uninitialized(1)) {
        assert(In.size() == 1);
        if constexpr (is_pair<typename Range::value_type>::value) {
            static_assert(std::is_pointer<typename Range::value_type::second_type>::value);
            parlay::assign_dispatch(pts[0], *(In[0].second), PointAssignTag());
        } else {
            parlay::assign_dispatch(pts[0], In[0], PointAssignTag());
        }
    }

    bool is_dummy;
    Points pts;
};

template<typename Point, typename AugType>
struct Interior : Node {
    Node* left;
    Node* right;
    AugType aug;
    Interior(Node* _left, Node* _right, AugType& _aug) :
        Node{false, _left->size + _right->size}, left(_left), right(_right), aug(_aug) {}
};

template<typename Point, typename AugType>
static Interior<Point, AugType>* allocInteriorNode(Node* L, Node* R, AugType& bit) {
    using Interior = Interior<Point, AugType>;
    Interior* o = parlay::type_allocator<Interior>::alloc();
    new (o) Interior(L, R, bit);
    return o;
}

// NOTE: Point: how data is stored in the tree
// NOTE: slice: input slice
template<typename Point, typename Range, typename LeafAllocTag, typename PointAssignTag>
static Leaf<Point, Range, PointAssignTag>* allocLeafNode(Range In, const auto alloc_size) {
    using leaf = Leaf<Point, Range, PointAssignTag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, alloc_size, LeafAllocTag());
    assert(o->is_dummy == false);
    return o;
}

template<typename Point, typename Range, typename LeafAllocTag, typename PointAssignTag>
static Leaf<Point, Range, PointAssignTag>* allocLeafNode(Range In) {
    using leaf = Leaf<Point, Range, PointAssignTag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, LeafAllocTag());
    assert(o->is_dummy == false);
    return o;
}

// template<typename Point, typename slice, typename assign_func>
// static leaf<Point, slice, assign_func>* alloc_dummy_leaf(
//     slice In,
//     assign_func const& ASSIGN_FUNC = [](Point& p, Point& q) { p = q; }) {
//   using leaf = leaf<Point, slice, assign_func>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf(In, true, ASSIGN_FUNC);
//   assert(o->is_dummy == true);
//   return o;
// }

// template<typename Point, typename slice, typename assign_function>
// static leaf<Point, slice, assign_function>* alloc_empty_leaf() {
//   using leaf = leaf<Point, slice, assign_function>;
//   leaf* o = parlay::type_allocator<leaf>::alloc();
//   new (o) leaf();
//   assert(o->size == 0 && o->pts.size() == 0);
//   return o;
// }

// template<typename Point>
// static void
// free_leaf(node* T) {
//   parlay::type_allocator<leaf<Point>>::retire(static_cast<leaf<Point>*>(T));
// }

template<typename Point, typename node_type>
static void free_node(Node* T) {
    // TODO: add static type check
    parlay::type_allocator<node_type>::retire(static_cast<node_type*>(T));
}

}  // namespace cpdd
