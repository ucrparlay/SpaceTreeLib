#ifndef PSTP_DEPENDENCE_TREE_NODE_H_
#define PSTP_DEPENDENCE_TREE_NODE_H_

#include <parlay/slice.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <tuple>
#include <type_traits>

#include "basic_point.h"
#include "comparator.h"
#include "parlay/utilities.h"

namespace pstp {

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
          typename PointAssignTag = parlay::move_assign_tag,
          bool kContainBox = false>
struct LeafNode : Node {
  using Points = parlay::sequence<Point>;
  using Box = std::pair<Point, Point>;  // TODO: use the version from Base tree

  static constexpr auto ContainBox() { return kContainBox; }

  constexpr auto InitBox() const {
    if constexpr (kContainBox) {
      return Box{};
    } else {
      return std::monostate{};
    }
  }

  // NOTE: default allocator
  LeafNode()
      : Node{true, static_cast<size_t>(0)}, is_dummy(false), split(InitBox()) {}

  // NOTE: alloc a leaf with default size
  LeafNode(Range In, AllocNormalLeafTag)
    requires(!kContainBox)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)) {
    assert(In.size() <= kDefaultWrap);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
  }

  // NOTE: alloc a leaf with default size and bounding box
  LeafNode(Range In, Box const& box, AllocNormalLeafTag)
    requires(kContainBox)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(kDefaultWrap)),
        split(box) {
    assert(In.size() <= kDefaultWrap);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
  }

  // NOTE: alloc a leaf with specific size
  LeafNode(Range In, size_t const alloc_size, AllocNormalLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(false),
        pts(Points::uninitialized(alloc_size)) {
    assert(In.size() <= alloc_size);
    std::ranges::for_each(In, [&, i = 0](auto&& x) mutable {
      parlay::assign_dispatch(pts[i++], x, PointAssignTag());
    });
  }

  // NOTE: alloc a dummy leaf
  LeafNode(Range In, AllocDummyLeafTag)
      : Node{true, static_cast<size_t>(In.size())},
        is_dummy(true),
        pts(Points::uninitialized(1)) {
    parlay::assign_dispatch(pts[0], In[0], PointAssignTag());
    if constexpr (kContainBox) {  // handling the boundingbox
      split = Box{pts[0], pts[0]};
    } else {
      (void)split;
    }
  }

  bool is_dummy;
  Points pts;
  std::conditional_t<kContainBox, Box, std::monostate> split;  // aka. MBR
};

// NOTE:: Alloc a leaf with input IN and given size
template <typename Range, typename Leaf>
static Leaf* AllocFixSizeLeafNode(Range In, size_t const alloc_size) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, alloc_size, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a leaf
template <typename Range, typename Leaf>
  requires(!Leaf::ContainBox())
static Leaf* AllocNormalLeafNode(Range In) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a leaf with bounding box
template <typename Range, typename Leaf, typename Box>
  requires(Leaf::ContainBox())
static Leaf* AllocNormalLeafNode(Range In, Box const& box) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf(In, box, AllocNormalLeafTag());
  assert(o->is_dummy == false);
  return o;
}

// NOTE: Alloc a empty Leaf
template <typename Range, typename Leaf>
  requires(!Leaf::ContainBox())
static Leaf* AllocEmptyLeafNode() {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf();
  return o;
}

// NOTE: Alloc a empty Leaf with bounding box
template <typename Range, typename Leaf, typename Box>
  requires(Leaf::ContainBox())
static Leaf* AllocEmptyLeafNode(Box const& box) {
  Leaf* o = parlay::type_allocator<Leaf>::alloc();
  new (o) Leaf();
  static_cast<Leaf*>(o)->split = box;
  return o;
}

// NOTE: Alloc a dummy Leaf (with bounding box)
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

  Node* left;
  Node* right;
  ST split;
  AT aug;
};

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

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_TREE_NODE_H_
