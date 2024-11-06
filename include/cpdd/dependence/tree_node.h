#pragma once
#include <parlay/slice.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <numeric>
#include <type_traits>

#include "basic_point.h"
#include "comparator.h"
#include "cpdd/base_tree.h"
#include "parlay/utilities.h"

namespace cpdd {

struct AllocNormalLeafTag {};
struct AllocDummyLeafTag {};
struct AllocEmptyLeafTag {};

struct Node {
  Node() : is_leaf{false}, size{0} {};
  Node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

  bool is_leaf;
  size_t size;
};

template <typename Point, typename Range, uint_fast8_t kDefaultWrap,
          typename PointAssignTag = parlay::move_assign_tag>
struct LeafNode : Node {
  using Points = parlay::sequence<Point>;

  // NOTE: default allocator
  LeafNode() : Node{true, static_cast<size_t>(0)}, is_dummy(false) {}

  // NOTE: alloc a normal leaf with size of DefaultWrap
  LeafNode(Range In, AllocNormalLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)) {
    assert(In.size() <= kDefaultWrap);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
    // for (int i = 0; i < In.size(); i++) {
    //     parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
    // }
  }

  // NOTE: alloc a normal leaf with specific size
  LeafNode(Range In, size_t const alloc_size, AllocNormalLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(alloc_size)) {
    assert(In.size() <= alloc_size);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
    // for (int i = 0; i < In.size(); i++) {
    //     parlay::assign_dispatch(pts[i], In[i], PointAssignTag());
    // }
  }

  // NOTE: alloc a dummy leaf
  LeafNode(Range In, AllocDummyLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(true),
        pts(Points::uninitialized(1)) {
    parlay::assign_dispatch(pts[0], In[0], PointAssignTag());
  }

  bool is_dummy;
  Points pts;
};

// NOTE:: Alloc a leaf with input IN and given size
template <typename Range, typename Leaf>
static Leaf* AllocFixSizeLeafNode(Range In, size_t const alloc_size) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, alloc_size, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a leaf with input IN and default leaf wrap
template <typename Range, typename Leaf>
static Leaf* AllocNormalLeafNode(Range In) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a empty Leaf
template <typename Range, typename Leaf>
static Leaf* AllocEmptyLeafNode() {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf();
  return o;
}

// NOTE: Alloc a dummy Leaf
template <typename Range, typename Leaf>
static Leaf* AllocDummyLeafNode(Range In) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, AllocDummyLeafTag());
  assert(o->is_dummy == true);
  return o;
}

template <typename Point, typename SplitType, typename AugType>
struct BinaryNode : Node {
  using PT = Point;
  using ST = SplitType;
  using AT = AugType;

  BinaryNode(Node* _left, Node* _right, const ST& _split, const AT& _aug)
      : Node{false, _left->size + _right->size},
        left(_left),
        right(_right),
        split(_split),
        aug(_aug) {}

  Node* left;
  Node* right;
  ST split;
  AT aug;
};

template <typename Point, uint_fast8_t kMD, typename SplitType,
          typename AugType>
struct MultiNode : Node {
  static constexpr uint_fast8_t const kRegions = 1 << kMD;

  using BucketType = uint_fast8_t;
  using Coord = typename Point::Coord;

  using Num = Num_Comparator<Coord>;
  using NodeArr = std::array<Node*, kRegions>;
  using ST = SplitType;
  using AT = AugType;

  MultiNode(NodeArr const& _tree_nodes, const ST& _split, const AT& _aug)
      : Node{false,
             std::accumulate(
                 _tree_nodes.begin(), _tree_nodes.end(), static_cast<size_t>(0),
                 [](size_t acc, Node* n) -> size_t { return acc + n->size; })},
        tree_nodes(_tree_nodes),
        split(_split),
        aug(_aug) {}

  // NOTE: Generate the hyperplane sequences for the node
  template <typename HyperPlaneSeq>
  inline void GenerateHyperPlaneSeq(HyperPlaneSeq& hyper_seq, auto idx,
                                    auto deep) {
    if (idx >= kRegions) {
      return;
    }
    hyper_seq[idx] = split[deep];
    GenerateHyperPlaneSeq(hyper_seq, idx * 2, deep + 1);
    GenerateHyperPlaneSeq(hyper_seq, idx * 2 + 1, deep + 1);
    return;
  }

  // NOTE: given an idx in the tree skeleton, modify its value its position in
  // left or right of the splitter
  inline BucketType SeievePoint(Point const& p, BucketType idx) {
    for (BucketType i = 0; i < kMD; ++i) {
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[i].second], split[i].first));
    }
    return idx;
  }

  inline size_t ReduceSums(BucketType const idx) {
    if (idx == 1) {
      return this->size;
    } else if (idx >= kRegions) {
      return tree_nodes[idx - kRegions]->size;
    } else {
      return ReduceSums(2 * idx) + ReduceSums(2 * idx + 1);
    }
  }

  template <typename BoxSeqSlice>
  void ComputeSubregionsRec(BoxSeqSlice box_seq, BucketType deep) {
    if (deep == kMD) {
      assert(box_seq.size() == 1);
      return;
    }

    BucketType n = box_seq.size();
    assert(n - n / 2 == n / 2);

    std::for_each_n(box_seq.begin(), n / 2, [&](auto& bx) {
      bx.second.pnt[split[deep].second] = split[deep].first;
    });
    std::for_each_n(box_seq.begin() + n / 2, n - n / 2, [&](auto& bx) {
      bx.first.pnt[split[deep].second] = split[deep].first;
    });

    ComputeSubregionsRec(box_seq.cut(0, n / 2), deep + 1);
    ComputeSubregionsRec(box_seq.cut(n / 2, n), deep + 1);
    return;
  }

  // NOTE: Given a box, produce new sub-boxes equivalent to modify the existing
  // boxes
  template <typename BoxSeq, typename Box>
  BoxSeq ComputeSubregions(Box const& box) {
    auto box_seq = BoxSeq(kRegions, box);
    ComputeSubregionsRec(parlay::make_slice(box_seq), 0);
    return std::move(box_seq);
  }

  // NOTE: Given a box and a bucket id, construct new box for that bucket
  template <typename Box>
  inline Box GetBoxByRegionId(BucketType id, Box const& box) {
    Box bx(box);
    assert(id >= 0 && id < kRegions);

    // PERF: cannot set i>=0 as it is unsigned int. idx 9 -> 101 -> RLR
    for (BucketType i = kMD; i > 0; --i) {
      auto& target = (id & (1 << (i - 1))) ? bx.first : bx.second;
      target.pnt[split[kMD - i].second] = split[kMD - i].first;
    }

    return std::move(bx);
  }

  // NOTE: Given a box and a bucket id, modify the box for that bucket
  template <typename Box>
  inline void ModifyBoxById(BucketType id, Box& box) {
    assert(id >= 0 && id <= kRegions);
    for (BucketType i = kMD; i > 0; --i) {
      auto& target = (id & (1 << (i - 1))) ? box.first : box.second;
      target.pnt[split[kMD - i].second] = split[kMD - i].first;
    }
    return;
  }

  NodeArr tree_nodes;
  ST split;
  AT aug;
};

template <typename Interior>
static Interior* AllocInteriorNode(Node* L, Node* R,
                                   typename Interior::ST const& split,
                                   typename Interior::AT const& aug) {
  Interior* o = parlay::type_allocator<Interior>::alloc();
  new (o) Interior(L, R, split, aug);
  return o;
}

template <typename Interior>
static Interior* AllocInteriorNode(typename Interior::NodeArr const& tree_nodes,
                                   typename Interior::ST const& split,
                                   typename Interior::AT const& aug) {
  Interior* o = parlay::type_allocator<Interior>::alloc();
  new (o) Interior(tree_nodes, split, aug);
  return o;
}

template <typename NodeType>
static void FreeNode(Node* T) {
  parlay::type_allocator<NodeType>::retire(static_cast<NodeType*>(T));
}

}  // namespace cpdd
