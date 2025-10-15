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
          typename AugType>
struct Node {
  using Coord = typename Point::Coord;
  using Num = Num_Comparator<Coord>;
  using ST = SplitType;
  using AT = AugType;

  static consteval auto GetRegions() { return 1 << kLevels; }

  static consteval auto GetLevels() { return kLevels; }

  Node() = default;

  Node(uint32_t _size, uint32_t _leaf_offset, ST const& _split, AT const& _aug,
       bool _is_leaf)
      : size(_size),
        leaf_offset(_leaf_offset),
        split(_split),
        aug(_aug),
        is_leaf(_is_leaf) {}

  // NOTE: up tp 4 billion points in 1.5 TB memory
  uint32_t size;         // Number of points in the subtree
  uint32_t leaf_offset;  // Offset to leaf points in leaf_seq
  ST split;              // Splitters for each dimension
  AT aug;                // Augmented data
  bool is_leaf;
};

}  // namespace array_view
}  // namespace psi

#endif  // PSI_DEPENDENCE_TREE_NODE_ARRAY_H_
