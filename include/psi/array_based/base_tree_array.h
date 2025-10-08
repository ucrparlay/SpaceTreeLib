#ifndef PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_
#define PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_

#include <cstdint>
#include <vector>
#include <limits>

#include "../dependence/comparator.h"
#include "../dependence/loggers.h"
#include "../dependence/search_container.h"
#include "../dependence/tree_node.h"

namespace psi {
namespace array_based {

//==============================================================================
// ARRAY-BASED BASE TREE CLASS
//
// A template-based spatial tree framework using array indices instead of
// pointers for better cache locality and reduced memory overhead.
//
// Key differences from pointer-based BaseTree:
// - Uses NodeIndex (uint32_t) instead of Node* pointers
// - Stores nodes in contiguous std::vector
// - Better cache locality for traversal operations
// - Lower memory overhead (4 bytes vs 8 bytes per reference)
//==============================================================================

template <typename Point, typename DerivedTree = void,
          uint_fast8_t kSkHeight = 6, uint_fast8_t kImbaRatio = 30>
class BaseTreeArray {
 public:
  //============================================================================
  // SECTION 1: TYPE DEFINITIONS & ALIASES
  //============================================================================

  // 1.1 Node Index Type (replaces Node*)
  using NodeIndex = uint32_t;
  static constexpr NodeIndex NULL_INDEX = std::numeric_limits<NodeIndex>::max();

  // 1.2 Point Types (same as pointer-based)
  using BasicPoint = typename Point::BP;
  using TemplatePoint = Point;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using DisType = typename Point::DisType;

  // 1.3 Dimension & Index Types (same as pointer-based)
  using DimsType = uint_fast8_t;
  using DepthType = int;
  using BucketType =
      std::conditional_t<(kSkHeight > 7), uint_fast16_t, uint_fast8_t>;
  using BallsType = uint_fast32_t;
  using IDType = uint_fast32_t;

  // 1.4 Container Types
  using Num = Num_Comparator<Coord>;
  using Slice = parlay::slice<Point*, Point*>;
  using ConstSlice = parlay::slice<Point const*, Point const*>;
  using Points = parlay::sequence<Point>;
  using ConstPoints = parlay::sequence<Point> const;
  using PointsIter = typename parlay::sequence<Point>::iterator;
  using BucketSeq = parlay::sequence<BucketType>;
  using BallSeq = parlay::sequence<BallsType>;

  // 1.5 Spatial Structures
  using HyperPlane = std::pair<Coord, DimsType>;
  using HyperPlaneSeq = parlay::sequence<HyperPlane>;
  using Box = std::pair<BasicPoint, BasicPoint>;
  using BoxSeq = parlay::sequence<Box>;

  // 1.6 Node-related Types (using indices instead of pointers)
  using NodeBoolean = std::pair<NodeIndex, bool>;
  using NodeBox = std::pair<NodeIndex, Box>;
  using NodeBoxSeq = parlay::sequence<NodeBox>;
  using NodeTag = std::pair<NodeIndex, uint_fast8_t>;
  using NodeTagSeq = parlay::sequence<NodeTag>;

  //============================================================================
  // SECTION 2: COMPILE-TIME CONSTANTS (same as pointer-based)
  //============================================================================

  static constexpr DimsType const kDim = std::tuple_size_v<Coords>;
  static constexpr BucketType const kBuildDepthOnce = kSkHeight;
  static constexpr BucketType const kPivotNum = (1 << kBuildDepthOnce) - 1;
  static constexpr BucketType const kBucketNum = 1 << kBuildDepthOnce;

  // Leaf configuration (same as pointer-based)
#ifdef LEAF1
  static constexpr uint_fast8_t const kLeaveWrap = 1;
#elif defined(LEAF2)
  static constexpr uint_fast8_t const kLeaveWrap = 2;
#elif defined(LEAF4)
  static constexpr uint_fast8_t const kLeaveWrap = 4;
#elif defined(LEAF16)
  static constexpr uint_fast8_t const kLeaveWrap = 16;
#elif defined(LEAF32)
  static constexpr uint_fast8_t const kLeaveWrap = 32;
#elif defined(LEAF64)
  static constexpr uint_fast8_t const kLeaveWrap = 64;
#elif defined(LEAF128)
  static constexpr uint_fast8_t const kLeaveWrap = 128;
#elif defined(LEAF512)
  static constexpr uint_fast16_t const kLeaveWrap = 512;
#elif defined(LEAF1024)
  static constexpr uint_fast16_t const kLeaveWrap = 1024;
#else
  static constexpr uint_fast8_t const kLeaveWrap = 32;  // Default
#endif

  static constexpr uint_fast8_t const kThinLeaveWrap = 24;
  static constexpr uint_fast8_t const kSlimLeaveWrap = 8;
  static constexpr uint_fast16_t const kSerialBuildCutoff = 1 << 10;
  static constexpr uint_fast8_t const kLog2Base = 10;
  static constexpr uint_fast16_t const kBlockSize = 1 << kLog2Base;
  static constexpr uint_fast8_t const kInbalanceRatio = kImbaRatio;

  //============================================================================
  // SECTION 3: VIRTUAL INTERFACE
  //============================================================================

  constexpr virtual void DeleteTree() = 0;

  virtual ~BaseTreeArray() = default;

  //============================================================================
  // SECTION 4: PUBLIC API (TO BE IMPLEMENTED BY DERIVED CLASSES)
  //============================================================================

  // Size and statistics
  size_t GetSize() const { return num_points_; }
  size_t GetNumNodes() const { return nodes_.size(); }
  size_t GetMemoryUsage() const {
    return nodes_.capacity() * sizeof(typename DerivedTree::Node) +
           leaf_points_.capacity() * sizeof(Point);
  }

  // Box operations (similar to pointer-based)
  static Box GetEmptyBox() {
    return Box(BasicPoint::GetMax(), BasicPoint::GetMin());
  }

  template <typename Range>
  static Box GetBox(Range const& In);

  static Box GetBox(Box const& a, Box const& b);

 protected:
  //============================================================================
  // SECTION 5: MEMBER VARIABLES
  //============================================================================

  NodeIndex root_ = NULL_INDEX;           // Root node index
  parlay::internal::timer timer;
  Box tree_box_;
  size_t num_points_ = 0;
  size_t num_nodes_ = 0;

  // Note: Actual node storage is in derived classes
  // Each derived class defines its own node structure
};

}  // namespace array_based
}  // namespace psi

// Include implementation files
#include "base_tree_array_impl/box_op.hpp"

#endif  // PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_
