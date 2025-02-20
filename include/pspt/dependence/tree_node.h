#ifndef PSPT_DEPENDENCE_TREE_NODE_H_
#define PSPT_DEPENDENCE_TREE_NODE_H_

#include <parlay/slice.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <variant>

#include "basic_point.h"
#include "comparator.h"
#include "parlay/utilities.h"

namespace pspt {

struct AllocNormalLeafTag {};
struct AllocDummyLeafTag {};
struct AllocEmptyLeafTag {};

struct Node {
  Node() : is_leaf{false}, size{0} {};
  Node(bool _is_leaf, size_t _size) : is_leaf{_is_leaf}, size{_size} {};

  // Adding a virtual destructor makes Node polymorphic
  virtual ~Node() = default;

  bool is_leaf;
  size_t size;
};

template <typename Point, typename Range, uint_fast8_t kDefaultWrap,
          typename AugType = std::monostate,
          typename PointAssignTag = parlay::move_assign_tag>
struct LeafNode : Node {
  using Points = parlay::sequence<Point>;
  using AT = AugType;

  // NOTE: default allocator
  LeafNode()
      : Node{true, static_cast<size_t>(0)}, is_dummy(false), aug(AugType()) {}

  LeafNode(AT const& _aug)
      : Node{true, static_cast<size_t>(0)}, is_dummy(false), aug(_aug) {}

  // NOTE: alloc a leaf with default size
  LeafNode(Range In, AugType const& _aug, AllocNormalLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)),
        aug(_aug) {
    assert(In.size() <= kDefaultWrap);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
  }

  // NOTE: alloc a leaf with specific size
  LeafNode(Range In, size_t const alloc_size, AugType const& _aug,
           AllocNormalLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(alloc_size)),
        aug(_aug) {
    assert(In.size() <= alloc_size);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
  }

  // NOTE: alloc a dummy leaf
  LeafNode(Range In, AugType const& _aug, AllocDummyLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(true),
        pts(Points::uninitialized(1)),
        aug(_aug) {
    parlay::assign_dispatch(pts[0], In[0], PointAssignTag());
  }

  bool CapacityFull() const { return this->size == kDefaultWrap; }

  Points const& GetPoints() const { return pts; }

  AugType const& GetAug() const { return aug; }

  AugType const& GetSplit() const { return aug; }

  bool is_dummy;
  Points pts;
  AugType aug;
};

// NOTE:: Alloc a leaf with input IN and given size
template <typename Range, typename Leaf>
static Leaf* AllocFixSizeLeafNode(
    Range In, size_t const alloc_size,
    typename Leaf::AT const& aug = typename Leaf::AT()) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, alloc_size, aug, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a leaf
template <typename Range, typename Leaf>
static Leaf* AllocNormalLeafNode(
    Range In, typename Leaf::AT const& aug = typename Leaf::AT()) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, aug, AllocNormalLeafTag());
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

// NOTE: Alloc a empty Leaf, but set the aug by hand
template <typename Range, typename Leaf>
static Leaf* AllocEmptyLeafNode(typename Leaf::AT const& aug) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(aug);
  return o;
}

// NOTE: Alloc a dummy leaf
template <typename Range, typename Leaf>
static Leaf* AllocDummyLeafNode(
    Range In, typename Leaf::AT const& aug = typename Leaf::AT()) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, aug, AllocDummyLeafTag());
  assert(o->is_dummy == true);
  assert(o->pts.size() == 1);
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

  // Adding a virtual destructor makes Node polymorphic
  virtual ~BinaryNode() override = default;

  // NOTE: test whether we can fetch @depth levels from @T
  template <typename Interior>
  static inline bool TestDepth(Node* T, int cur_depth, int depth) {
    if (cur_depth == depth) {
      return true;
    }
    if (T->is_leaf) {
      return false;
    }
    auto TI = static_cast<Interior*>(T);
    return TestDepth<Interior>(TI->left, cur_depth + 1, depth) &&
           TestDepth<Interior>(TI->right, cur_depth + 1, depth);
  }

  ST const& GetSplit() const { return split; }

  Node* left;
  Node* right;
  ST split;
  AT aug;
};

// NOTE: multi-way node whose children # is fixed
template <typename Point, uint_fast8_t kMD, typename SplitType,
          typename AugType>
struct MultiNode : Node {
  using BucketType = uint_fast8_t;
  using Coord = typename Point::Coord;
  using Num = Num_Comparator<Coord>;
  using ST = SplitType;
  using AT = AugType;

  static consteval auto GetRegions() { return 1 << kMD; }

  static consteval auto GetLevels() { return kMD; }

  static consteval auto GetSplitNums() { return kMD; }

  // NOTE: whether we use same splitter for same level
  static consteval auto EqualSplit()
    requires std::same_as<ST, std::array<typename ST::value_type, kMD>>
  {
    return true;  // every dimension one splitter
  }

  static consteval auto EqualSplit()
    requires std::same_as<ST, std::array<typename ST::value_type, 1 << kMD>>
  {
    return false;  // the spliiter number mathes inter nodes num
  }

  using NodeArr = std::array<Node*, GetRegions()>;

  MultiNode(NodeArr const& _tree_nodes, const ST& _split, const AT& _aug)
      : Node{false,
             std::accumulate(
                 _tree_nodes.begin(), _tree_nodes.end(), static_cast<size_t>(0),
                 [](size_t acc, Node* n) -> size_t { return acc + n->size; })},
        tree_nodes(_tree_nodes),
        split(_split),
        aug(_aug) {}

  // Adding a virtual destructor makes Node polymorphic
  // TODO: check whether it is possible to remove it then KNN-mix
  virtual ~MultiNode() override = default;

  inline size_t MergeSize(BucketType const idx) {
    if (idx == 1) {
      return this->size;
    } else if (idx >= GetRegions()) {
      return tree_nodes[idx - GetRegions()]->size;
    } else {
      return MergeSize(2 * idx) + MergeSize(2 * idx + 1);
    }
  }

  ST const& GetSplit() const { return split; }

  NodeArr tree_nodes;
  ST split;
  AT aug;
};

// NOTE: dynamic multi-way node whose children # is dynamic
template <typename Point, typename SplitType, typename AugType>
struct DynamicNode : Node {
  using BucketType = uint_fast8_t;
  using Coord = typename Point::Coord;
  using Num = Num_Comparator<Coord>;
  using ST = SplitType;
  using AT = AugType;
  using NodeArr = parlay::sequence<Node*>;

  DynamicNode()
      : Node{false, static_cast<size_t>(0)},
        tree_nodes(NodeArr()),
        split(ST()),
        aug(AT()) {}

  DynamicNode(NodeArr const& _tree_nodes, const ST& _split, const AT& _aug)
      : Node{false,
             std::accumulate(
                 _tree_nodes.begin(), _tree_nodes.end(), static_cast<size_t>(0),
                 [](size_t acc, Node* n) -> size_t { return acc + n->size; })},
        tree_nodes(_tree_nodes),
        split(_split),
        aug(_aug) {}

  // Adding a virtual destructor makes Node polymorphic
  // TODO: check whether it is possible to remove it then KNN-mix
  virtual ~DynamicNode() override = default;

  ST const& GetSplit() const { return split; }

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

// TODO: maybe replace the NodeArr, ST, and AT with template parameters
template <typename Interior>
static Interior* AllocInteriorNode(typename Interior::NodeArr const& tree_nodes,
                                   typename Interior::ST const& split,
                                   typename Interior::AT const& aug) {
  Interior* o = parlay::type_allocator<Interior>::alloc();
  new (o) Interior(tree_nodes, split, aug);
  return o;
}

template <typename Interior>
static Interior* AllocInteriorNode() {
  Interior* o = parlay::type_allocator<Interior>::alloc();
  new (o) Interior();
  return o;
}

template <typename NodeType>
static void FreeNode(Node* T) {
  parlay::type_allocator<NodeType>::retire(static_cast<NodeType*>(T));
}

}  // namespace pspt

#endif  // PSPT_DEPENDENCE_TREE_NODE_H_
