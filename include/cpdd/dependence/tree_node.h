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
         typename PointAssignTag = parlay::move_assign_tag>
struct LeafNode : Node {
    using Points = parlay::sequence<Point>;

    // NOTE: default allocator
    LeafNode() : Node{true, static_cast<size_t>(0)}, is_dummy(false) {}

    // NOTE: alloc a normal leaf with size of DefaultWrap
    LeafNode(Range In, AllocNormalLeafTag) :
        Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)) {
        assert(In.size() <= kDefaultWrap);
        for (int i = 0; i < In.size(); i++) {
            parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
        }
    }

    // NOTE: alloc a normal leaf with specific size
    LeafNode(Range In, const size_t alloc_size, AllocNormalLeafTag) :
        Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(alloc_size)) {
        assert(In.size() <= alloc_size);
        for (int i = 0; i < In.size(); i++) {
            parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
        }
    }

    // NOTE: alloc a dummy leaf
    LeafNode(Range In, AllocDummyLeafTag) :
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
template<typename Range, typename Leaf>
static Leaf* AllocFixSizeLeafNode(Range In, const size_t alloc_size) {
    Leaf* o = parlay::type_allocator<Leaf>::alloc();
    new (o) Leaf(In, alloc_size, AllocNormalLeafTag());
    assert(o->is_dummy == false);
    return o;
}

// NOTE: Alloc a leaf with input IN and default leaf wrap
template<typename Range, typename Leaf>
static Leaf* AllocNormalLeafNode(Range In) {
    Leaf* o = parlay::type_allocator<Leaf>::alloc();
    new (o) Leaf(In, AllocNormalLeafTag());
    assert(o->is_dummy == false);
    return o;
}

// NOTE: Alloc a empty Leaf
template<typename Range, typename Leaf>
static Leaf* AllocEmptyLeafNode() {
    Leaf* o = parlay::type_allocator<Leaf>::alloc();
    new (o) Leaf();
    return o;
}

// NOTE: Alloc a dummy Leaf
template<typename Range, typename Leaf>
static Leaf* AllocDummyLeafNode(Range In) {
    Leaf* o = parlay::type_allocator<Leaf>::alloc();
    new (o) Leaf(In, AllocDummyLeafTag());
    assert(o->is_dummy == false);
    return o;
}

template<typename Point, typename SplitType, typename AugType>
struct InteriorNode : Node {
    using ST = SplitType;
    using AT = AugType;

    InteriorNode(Node* _left, Node* _right, const ST& _split, const AT& _aug) :
        Node{false, _left->size + _right->size},
        left(_left),
        right(_right),
        split(_split),
        aug(_aug) {}

    Node* left;
    Node* right;
    ST split;
    AT aug;
};

template<typename Interior>
static Interior* AllocInteriorNode(Node* L, Node* R,
                                   const typename Interior::ST& split,
                                   const typename Interior::AT& aug) {
    Interior* o = parlay::type_allocator<Interior>::alloc();
    new (o) Interior(L, R, split, aug);
    return o;
}

template<typename NodeType>
static void FreeNode(Node* T) {
    // TODO: add static type check
    parlay::type_allocator<NodeType>::retire(static_cast<NodeType*>(T));
}

}  // namespace cpdd
