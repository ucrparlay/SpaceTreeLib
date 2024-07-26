#pragma once
#include <parlay/slice.h>
#include <cstdint>
#include "basic_point.h"
#include "parlay/utilities.h"
#include "utility.h"

namespace cpdd {

struct Node {
    Node() : is_leaf{false}, size{0} {};
    Node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

    bool is_leaf;
    size_t size;
};

template<typename Point, typename Range, uint_fast8_t kDefaultWrap,
         typename PointAssignTag>
struct Leaf : Node {
    using Points = parlay::sequence<Point>;

    // NOTE: default allocator
    Leaf() : Node{true, static_cast<size_t>(0)}, is_dummy(false) {}

    // NOTE: alloc a normal leaf with size of DefaultWrap
    Leaf(Range In, AllocNormalLeafTag) :
        Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)) {
        assert(In.size() <= kDefaultWrap);
        for (int i = 0; i < In.size(); i++) {
            parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
        }
    }

    // NOTE: alloc a normal leaf with specific size
    Leaf(Range In, const auto alloc_size, AllocNormalLeafTag) :
        Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(alloc_size)) {
        assert(In.size() <= alloc_size);
        for (int i = 0; i < In.size(); i++) {
            parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
        }
    }

    // NOTE: alloc a dummy leaf
    Leaf(Range In, AllocDummyLeafTag) :
        Node{true, static_cast<size_t>(In.size())},
        is_dummy(true),
        pts(Points::uninitialized(1)) {
        assert(In.size() == 1);
        parlay::assign_dispatch(pts[0], In[0], PointAssignTag());
    }

    bool is_dummy;
    Points pts;
};

// NOTE:: Alloc a leaf with input IN and given size
template<typename Point, typename Range, typename LeafAllocTag,
         uint_fast8_t kDefaultWrap,
         typename PointAssignTag = parlay::move_assign_tag>
static Leaf<Point, Range, kDefaultWrap, PointAssignTag>* AllocLeafNode(
    Range In, const auto alloc_size) {
    static_assert(std::is_integral<decltype(alloc_size)>::value);

    using leaf = Leaf<Point, Range, kDefaultWrap, PointAssignTag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, alloc_size, LeafAllocTag());
    assert(o->is_dummy == false);
    return o;
}

// NOTE: Alloc a leaf with input IN and default leaf wrap
template<typename Point, typename Range, typename LeafAllocTag,
         uint_fast8_t kDefaultWrap,
         typename PointAssignTag = parlay::move_assign_tag>
static Leaf<Point, Range, kDefaultWrap, PointAssignTag>* AllocLeafNode(
    Range In) {
    using leaf = Leaf<Point, Range, kDefaultWrap, PointAssignTag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf(In, LeafAllocTag());
    assert(o->is_dummy == false);
    return o;
}

// NOTE: Alloc a empty leaf
template<typename Point, typename Range, typename LeafAllocTag,
         uint_fast8_t kDefaultWrap,
         typename PointAssignTag = parlay::move_assign_tag>
static Leaf<Point, Range, kDefaultWrap, PointAssignTag>* AllocLeafNode() {
    using leaf = Leaf<Point, Range, kDefaultWrap, PointAssignTag>;
    leaf* o = parlay::type_allocator<leaf>::alloc();
    new (o) leaf();
    return o;
}

template<typename Point, typename SplitType, typename AugType>
struct Interior : Node {
    Node* left;
    Node* right;
    SplitType split;
    AugType aug;
    Interior(Node* _left, Node* _right, const SplitType& _split,
             const AugType& _aug) :
        Node{false, _left->size + _right->size},
        left(_left),
        right(_right),
        split(_split),
        aug(_aug) {}
};

template<typename Point, typename SplitType, typename AugType>
static Interior<Point, SplitType, AugType>* AllocInteriorNode(
    Node* L, Node* R, const SplitType& split, const AugType& aug) {
    using Interior = Interior<Point, SplitType, AugType>;
    Interior* o = parlay::type_allocator<Interior>::alloc();
    new (o) Interior(L, R, split, aug);
    return o;
}

template<typename Point, typename NodeType>
static void FreeNode(Node* T) {
    // TODO: add static type check
    parlay::type_allocator<NodeType>::retire(static_cast<NodeType*>(T));
}

}  // namespace cpdd
