#ifndef PSI_DEPENDENCE_TREE_NODE_ARRAY_H_
#define PSI_DEPENDENCE_TREE_NODE_ARRAY_H_

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
#include "dependence/concepts.h"
#include "parlay/utilities.h"

namespace psi {
namespace array_view {
//============================================================================
// ARRAY-BASED NODE STRUCTURES
//============================================================================

template <typename Point, uint_fast8_t kLevels, typename SplitType,
          typename AugType, uint_fast8_t kLeafWrap, typename PointAssignTagType>
struct ArrayNode {
  using Coord = typename Point::Coord;
  using Num = Num_Comparator<Coord>;
  using ST = SplitType;
  using AT = AugType;
  using NodeIndex = uint32_t;  // TODO: use the one in the type_traits
  using PointAssignTag = PointAssignTagType;

  // TODO: Add how to access its children

  static consteval auto GetRegions() { return 1 << kLevels; }

  static consteval auto GetLevels() { return kLevels; }

  static consteval auto GetLeafWrap() { return kLeafWrap; }

  ArrayNode()
      : size(0),
        leaf_offset(0),
        split(),
        aug(),
        is_leaf(false),
        is_dummy(false) {}

  ArrayNode(NodeIndex _size, NodeIndex _leaf_offset, ST const& _split,
            AT const& _aug, bool _is_leaf, bool _is_dummy)
      : size(_size),
        leaf_offset(_leaf_offset),
        split(_split),
        aug(_aug),
        is_leaf(_is_leaf),
        is_dummy(_is_dummy) {}

  void Update(NodeIndex _size, NodeIndex _leaf_offset, ST const& _split,
              AT const& _aug, bool _is_leaf) {}

  // NOTE: up tp 4 billion 2-d points in 1.5 TB memory
  NodeIndex size;         // Number of points in the subtree
  NodeIndex leaf_offset;  // Offset to leaf points in leaf_seq
  ST split;               // Splitters for each dimension
  AT aug;                 // Augmented data
  bool is_leaf;           // Whether this points to a leaf node
  bool is_dummy = false;  // Whether this is a dummy leaf
};

template <typename InnerSeq, typename NodeIndex>
void AnnoteEmptyLeaf(InnerSeq& inner_tree_seq, NodeIndex inner_offset,
                     NodeIndex leaf_offset) {
  using Node = typename InnerSeq::value_type;

  inner_tree_seq[inner_offset].size = 0;
  inner_tree_seq[inner_offset].leaf_offset = leaf_offset;
  inner_tree_seq[inner_offset].is_leaf = true;
  // ignore aug and split for leaf
  assert(!HasBox<typename Node::AT> ||
         inner_tree_seq[inner_offset].GetBox() == Node::AT::Geo::GetEmptyBox());
  return;
}

template <typename InnerSeq, typename LeafSeq, typename InputSlice,
          typename NodeIndex>
void AnnoteLeaf(InnerSeq& inner_tree_seq, LeafSeq& leaf_seq, InputSlice In,
                NodeIndex inner_offset, NodeIndex leaf_offset) {
  using Node = typename InnerSeq::value_type;

  inner_tree_seq[inner_offset].size = In.size();
  inner_tree_seq[inner_offset].leaf_offset = leaf_offset;
  inner_tree_seq[inner_offset].is_leaf = true;
  inner_tree_seq[inner_offset].is_dummy = false;
  inner_tree_seq[inner_offset].aug = Node::AT::Create(In);
  // ignore split for leaf

  for (int i = 0; i < In.size(); i++) {
    parlay::assign_dispatch(leaf_seq[leaf_offset + i], In[i],
                            Node::PointAssignTag());
  }
  return;
}

template <typename InnerSeq, typename LeafSeq, typename InputSlice,
          typename NodeIndex>
void AnnoteDummyLeaf(InnerSeq& inner_tree_seq, LeafSeq& leaf_seq, InputSlice In,
                     NodeIndex inner_offset, NodeIndex leaf_offset) {
  using Node = typename InnerSeq::value_type;

  assert(In.size > 0);

  inner_tree_seq[inner_offset].size = In.size();
  inner_tree_seq[inner_offset].leaf_offset = leaf_offset;
  inner_tree_seq[inner_offset].is_leaf = true;
  inner_tree_seq[inner_offset].is_dummy = true;
  inner_tree_seq[inner_offset].aug = Node::AT::Create(In.cut(0, 1));

  // ignore aug and split for leaf
  parlay::assign_dispatch(leaf_seq[leaf_offset], In[0], Node::PointAssignTag());
  return;
}

template <typename InnerSeq, typename NodeIndex, typename SplitType,
          typename AugType>
void AnnoteInterior(InnerSeq& inner_tree_seq, NodeIndex inner_offset,
                    NodeIndex leaf_offset, NodeIndex size,
                    SplitType const split) {
  using Node = typename InnerSeq::value_type;

  inner_tree_seq[inner_offset].size = size;
  inner_tree_seq[inner_offset].leaf_offset = leaf_offset;
  inner_tree_seq[inner_offset].is_leaf = false;
  inner_tree_seq[inner_offset].is_dummy = false;
  inner_tree_seq[inner_offset].split = split;
  inner_tree_seq[inner_offset].aug =
      Node::AT::Create(inner_tree_seq, inner_offset);

  return;
}

}  // namespace array_view
}  // namespace psi

#endif  // PSI_DEPENDENCE_TREE_NODE_ARRAY_H_
