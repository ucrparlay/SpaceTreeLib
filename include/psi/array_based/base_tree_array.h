#ifndef PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_
#define PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_

#include <cstdint>
#include <limits>
#include <vector>

#include "dependence/comparator.h"
#include "dependence/geo_base.h"
#include "dependence/loggers.h"
#include "dependence/search_container.h"
#include "dependence/tree_node_array.h"

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

template <typename TypeTrait, typename DerivedTree = void>
class BaseTreeArray {
 public:
  //============================================================================
  // SECTION 1: TYPE DEFINITIONS & ALIASES
  //============================================================================

  // 1.1 Node Index Type (replaces Node*)
  using NodeIndex = uint32_t;
  static constexpr NodeIndex NULL_INDEX = std::numeric_limits<NodeIndex>::max();

  // 1.1 Point Types
  using BasicPoint = typename TypeTrait::BasicPoint;
  using Coord = typename TypeTrait::Coord;
  using Coords = typename TypeTrait::Coords;
  using DisType = typename TypeTrait::DisType;
  using Point = typename TypeTrait::Point;
  using TemplatePoint = typename TypeTrait::TemplatePoint;

  // 1.2 Dimension & Index Types
  using BallsType = typename TypeTrait::BallsType;
  using BucketType = typename TypeTrait::BucketType;
  using DepthType = typename TypeTrait::DepthType;
  using DimsType = typename TypeTrait::DimsType;
  using IDType = typename TypeTrait::IDType;

  // 1.3 Container Types
  using BallSeq = typename TypeTrait::BallSeq;
  using BucketSeq = typename TypeTrait::BucketSeq;
  using ConstPoints = typename TypeTrait::ConstPoints;
  using ConstSlice = typename TypeTrait::ConstSlice;
  using Num = typename TypeTrait::Num;
  using Points = typename TypeTrait::Points;
  using PointsIter = typename TypeTrait::PointsIter;
  using Slice = typename TypeTrait::Slice;

  // 1.4 Spatial Structures
  using Box = typename TypeTrait::Box;
  using BoxSeq = typename TypeTrait::BoxSeq;
  using HyperPlane = typename TypeTrait::HyperPlane;
  using HyperPlaneSeq = typename TypeTrait::HyperPlaneSeq;

  // 1.5 Node & Tree Structures
  using NodeBoolean = typename TypeTrait::NodeBoolean;
  using NodeBox = typename TypeTrait::NodeBox;
  using NodeBoxSeq = typename TypeTrait::NodeBoxSeq;
  using NodeTag = typename TypeTrait::NodeTag;
  using NodeTagSeq = typename TypeTrait::NodeTagSeq;

  // 1.6 Geometric Base
  using Geo = GeoBase<TypeTrait>;
  using NormalCircle = typename Geo::NormalCircle;
  using CoverCircle = typename Geo::CoverCircle;
  using BoxCut = typename Geo::BoxCut;

  //============================================================================
  // SECTION 2: COMPILE-TIME CONSTANTS (same as pointer-based)
  //============================================================================

  static constexpr DimsType const kDim = std::tuple_size_v<Coords>;
  static constexpr BucketType const kBuildDepthOnce = TypeTrait::kSkHeight;
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
  static constexpr uint_fast8_t const kInbalanceRatio = TypeTrait::kImbaRatio;

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
  auto GetRoot() const { return root_; }

  // Box operations (similar to pointer-based)
  static Box GetEmptyBox() {
    return Box(BasicPoint::GetMax(), BasicPoint::GetMin());
  }

  consteval static auto GetBuildDepthOnce() {
    return static_cast<int>(kBuildDepthOnce);
  }

  static inline size_t GetImbalanceRatio() {
    return static_cast<size_t>(kInbalanceRatio);
  }

 protected:
  //============================================================================
  // SECTION 5: MEMBER VARIABLES
  //============================================================================

  NodeIndex root_ = NULL_INDEX;  // Root node index
  parlay::internal::timer timer;
  Box tree_box_;
  size_t num_points_ = 0;
  size_t num_nodes_ = 0;

  // Note: Actual node storage is in derived classes
  // Each derived class defines its own node structure
};

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_BASE_TREE_ARRAY_H_
