#ifndef PSI_BASE_TREE_H_
#define PSI_BASE_TREE_H_

#include <sys/types.h>

#include <cstdint>
#include <cstdio>
#include <type_traits>

#include "dependence/comparator.h"
#include "dependence/loggers.h"
#include "dependence/search_container.h"
#include "dependence/tree_node.h"

namespace psi {

//==============================================================================
// BASE TREE CLASS
//
// A template-based spatial tree framework using CRTP (Curiously Recurring
// Template Pattern) for zero-cost abstraction. Supports orthogonal, k-d,
// and cover tree implementations.
//
// Implementation is split across base_tree_impl/ for maintainability.
// See REFACTORING_NOTES.md for design rationale.
//==============================================================================

template <typename Point, typename DerivedTree = void,
          uint_fast8_t kSkHeight = 6, uint_fast8_t kImbaRatio = 30>
class BaseTree {
 public:
  //============================================================================
  // SECTION 1: TYPE DEFINITIONS & ALIASES
  // Core type system used throughout the tree implementation
  //============================================================================

  // 1.1 Point Types
  using BasicPoint = Point::BP;
  using TemplatePoint = Point;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using DisType = typename Point::DisType;

  // 1.2 Dimension & Index Types
  using DimsType = uint_fast8_t;
  using DepthType = int;
  using BucketType =
      std::conditional_t<(kSkHeight > 7), uint_fast16_t, uint_fast8_t>;
  using BallsType = uint_fast32_t;
  using IDType = uint_fast32_t;

  // 1.3 Container Types
  using Num = Num_Comparator<Coord>;
  using Slice = parlay::slice<Point*, Point*>;
  using ConstSlice = parlay::slice<Point const*, Point const*>;
  using Points = parlay::sequence<Point>;
  using ConstPoints = parlay::sequence<Point> const;
  using PointsIter = typename parlay::sequence<Point>::iterator;
  using BucketSeq = parlay::sequence<BucketType>;
  using BallSeq = parlay::sequence<BallsType>;

  // 1.4 Spatial Structures
  using HyperPlane = std::pair<Coord, DimsType>;
  using HyperPlaneSeq = parlay::sequence<HyperPlane>;
  using Box = std::pair<BasicPoint, BasicPoint>;
  using BoxSeq = parlay::sequence<Box>;

  // 1.5 Node Types
  using NodeBoolean = std::pair<Node*, bool>;
  using NodeBox = std::pair<Node*, Box>;
  using NodeBoxSeq = parlay::sequence<NodeBox>;
  using NodeTag = std::pair<Node*, uint_fast8_t>;
  using NodeTagSeq = parlay::sequence<NodeTag>;

  //============================================================================
  // SECTION 2: COMPILE-TIME CONSTANTS
  // Configuration parameters determined at compile time
  //============================================================================

  static constexpr DimsType const kDim = std::tuple_size_v<Coords>;
  static constexpr BucketType const kBuildDepthOnce = kSkHeight;
  static constexpr BucketType const kPivotNum = (1 << kBuildDepthOnce) - 1;
  static constexpr BucketType const kBucketNum = 1 << kBuildDepthOnce;

  // Leaf configuration (controlled by compile-time flags)
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
#endif  // LEAF1

  static constexpr uint_fast8_t const kThinLeaveWrap = 24;
  static constexpr uint_fast8_t const kSlimLeaveWrap = 8;
  static constexpr uint_fast16_t const kSerialBuildCutoff = 1 << 10;
  static constexpr uint_fast8_t const kLog2Base = 10;
  static constexpr uint_fast16_t const kBlockSize = 1 << kLog2Base;
  static constexpr uint_fast8_t const kInbalanceRatio = kImbaRatio;

  //============================================================================
  // SECTION 3: NESTED TYPES
  // Forward declarations and nested class definitions
  //============================================================================

  template <typename Leaf, typename Interior>
  struct InnerTree;

  struct BoxCut;

  struct NormalCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return radius; }

    Coord GetRadiusSquare() const { return radius * radius; }

    static Coord ComputeRadius(double const& r) {
      return static_cast<Coord>(r);
    }

    friend std::ostream& operator<<(std::ostream& o, NormalCircle const& cl) {
      o << "{ " << cl.center << ", " << cl.radius << "}";
      return o;
    }

    Point center;
    Coord radius;
  };

  struct CoverCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return static_cast<Coord>(1 << level); }

    Coord GetRadiusSquare() const { return GetRadius() * GetRadius(); }

    static DepthType ComputeRadius(auto const& r) {
      return Num_Comparator<decltype(r)>::IsZero(r)
                 ? -1
                 : static_cast<DepthType>(std::ceil(std::log2(r)));
    }

    bool operator==(CoverCircle const& cl) const {
      return center == cl.center && level == cl.level;
    }

    friend std::ostream& operator<<(std::ostream& o, CoverCircle const& cl) {
      o << "{ " << cl.center << ", " << cl.level << "}";
      return o;
    }

    Point center;
    DepthType level;
  };

  //============================================================================
  // SECTION 4: PUBLIC API - CONSTRUCTION & LIFECYCLE
  // Implementation: base_tree_impl/delete_tree.hpp
  //============================================================================

  constexpr virtual void DeleteTree() = 0;

  template <typename Leaf, typename Interior>
  void DeleteTreeWrapper();

  template <typename Leaf, IsBinaryNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  template <typename Leaf, IsMultiNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  template <typename Leaf, IsDynamicNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  //============================================================================
  // SECTION 5: PUBLIC API - QUERY OPERATIONS
  // Implementation: base_tree_impl/knn_query.hpp
  //                 base_tree_impl/range_query.hpp
  //============================================================================

  // 5.1 K-Nearest Neighbor Search
  template <typename Leaf, typename Range>
  static void KNNLeaf(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq);

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void KNNBinary(Node* T, Point const& q,
                        kBoundedQueue<Point, Range>& bq, Box const& node_box,
                        KNNLogger& logger)
    requires std::same_as<typename Interior::ST, HyperPlane>;

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void KNNBinaryBox(Node* T, Point const& q,
                           kBoundedQueue<Point, Range>& bq, KNNLogger& logger);

  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMultiExpand(Node* T, Point const& q, DimsType dim,
                             BucketType idx, kBoundedQueue<Point, Range>& bq,
                             Box const& node_box, KNNLogger& logger);

  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMultiExpandBox(Node* T, Point const& q, DimsType dim,
                                BucketType idx, kBoundedQueue<Point, Range>& bq,
                                KNNLogger& logger);

  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMulti(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                       KNNLogger& logger);

  template <typename Leaf, IsBinaryNode BN, IsMultiNode MN, typename Range>
  static void KNNMix(Node* T, Point const& q, DimsType dim, BucketType idx,
                     kBoundedQueue<Point, Range>& bq, Box const& node_box,
                     KNNLogger& logger);

  template <typename Leaf, IsDynamicNode Interior, typename Range>
  static void KNNCover(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                       KNNLogger& logger);

  // 5.2 Range Query & Count
  template <typename Leaf>
  static size_t RangeCountRectangleLeaf(Node* T, Box const& query_box);

  template <typename Leaf, IsBinaryNode Interior>
  static size_t RangeCountRectangle(Node* T, Box const& query_box,
                                    Box const& node_box,
                                    RangeQueryLogger& logger);

  template <typename Leaf, IsMultiNode Interior>
  static size_t RangeCountRectangle(Node* T, Box const& query_box,
                                    Box const& node_box, DimsType dim,
                                    BucketType idx, RangeQueryLogger& logger);

  template <typename Leaf, typename Range>
  static void RangeQueryLeaf(Node* T, Range Out, size_t& s,
                             Box const& query_box);

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                        Box const& query_box,
                                        Box const& node_box,
                                        RangeQueryLogger& logger);

  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                        Box const& query_box,
                                        Box const& node_box, DimsType dim,
                                        BucketType idx,
                                        RangeQueryLogger& logger);

  //============================================================================
  // SECTION 6: PUBLIC API - TREE MODIFICATION OPERATIONS
  // Implementation: base_tree_impl/tree_op/rebuild.hpp
  //                 base_tree_impl/tree_op/leaf_op.hpp
  //============================================================================

  // 6.1 Insert Operations
  template <typename Leaf>
  static Node* InsertPoints2Leaf(Node* T, Slice In);

  template <typename Leaf, typename Interior, typename PrepareFunc,
            typename... Args>
  Node* RebuildWithInsert(Node* T, PrepareFunc prepare_func, Slice In,
                          Args&&... args);

  // 6.2 Delete Operations
  template <typename Leaf, typename RT>
  static RT DeletePoints4Leaf(Node* T, Slice In);

  template <typename Leaf, typename RT>
  static RT DiffPoints4Leaf(Node* T, Slice In);

  // 6.3 Rebuild Operations
  template <typename Leaf, typename Interior, bool granularity = true>
  static void PrepareRebuild(Node* T, Slice In, Points& wx, Points& wo);

  template <typename Leaf, typename Interior, bool granularity = true>
  static void PrepareRebuild(Node* T, Points& wx, Points& wo);

  template <typename Leaf, typename Interior, bool granularity,
            typename... Args>
  Node* RebuildSingleTree(Node* T, Args&&... args);

  template <typename Leaf, IsBinaryNode Interior, bool granularity,
            typename PrepareFunc, typename... Args>
  Node* RebuildTreeRecursive(Node* T, PrepareFunc&& prepare_func,
                             bool const allow_rebuild, Args&&... args);

  template <typename Leaf, IsMultiNode Interior, bool granularity,
            typename PrepareFunc, typename... Args>
  Node* RebuildTreeRecursive(Node* T, PrepareFunc&& prepare_func,
                             bool const allow_rebuild, Args&&... args);

  //============================================================================
  // SECTION 7: PUBLIC API - TREE INSPECTION & PROPERTIES
  // Implementation: base_tree_impl/validation.hpp
  //============================================================================

  void SetRoot(Node* root) { this->root_ = root; }

  Node* GetRoot() { return this->root_; }

  size_t GetSize() { return this->root_ ? this->root_->size : 0; }

  Box GetRootBox() { return this->tree_box_; }

  consteval static auto GetBuildDepthOnce() {
    return static_cast<int>(kBuildDepthOnce);
  }

  template <typename Leaf, typename Interior>
  size_t GetTreeHeight();

  template <typename Leaf, typename Interior>
  size_t GetMaxTreeDepth(Node* T, size_t deep);

  template <typename Leaf, typename Interior>
  double GetAveTreeHeight();

  template <typename Leaf, typename Interior>
  size_t CountTreeNodesNum(Node* T);

  template <typename Leaf, typename Interior>
  void CountTreeHeights(Node* T, size_t deep, size_t& idx,
                        parlay::sequence<size_t>& heights);

  //============================================================================
  // SECTION 8: PUBLIC API - VALIDATION & DEBUGGING
  // Implementation: base_tree_impl/validation.hpp
  //============================================================================

  template <typename Leaf, typename Interior>
  Box CheckBox(Node* T, Box const& box);

  template <typename Leaf, typename Interior>
  Points CheckCover(Node* T,
                    typename Interior::CircleType const& level_cover_circle);

  template <typename Leaf, typename Interior>
  static size_t CheckSize(Node* T);

  template <typename Leaf, typename Interior>
  void CheckTreeSameSequential(Node* T, int dim);

  template <typename Leaf, typename Interior, typename SplitRule>
  void Validate();

  //============================================================================
  // SECTION 9: GEOMETRY UTILITIES - BOX OPERATIONS
  // Implementation: base_tree_impl/box_op.hpp
  //============================================================================

  static inline Coord GetBoxMid(DimsType const d, Box const& bx);
  static inline bool LegalBox(Box const& bx);
  static inline bool WithinBox(Box const& a, Box const& b);
  static inline bool SameBox(Box const& a, Box const& b);
  static inline bool WithinBox(Point const& p, Box const& bx);
  static inline bool BoxIntersectBox(Box const& a, Box const& b);
  static inline bool IsBoxLineInDimension(Box const& box, DimsType d);
  static inline bool VerticalLineSplitBox(Coord const& l, Box const& box,
                                          DimsType d);
  static inline bool VerticalLineOnBoxLeftEdge(Coord const& l, Box const& box,
                                               DimsType d);
  static inline bool VerticalLineOnBoxRightEdge(Coord const& l, Box const& box,
                                                DimsType d);
  static inline bool VerticalLineOnBoxEdge(Coord const& l, Box const& box,
                                           DimsType d);
  static inline bool VerticalLineIntersectBox(Coord const& l, Box const& box,
                                              DimsType d);
  static inline bool VerticalLineIntersectBoxExclude(Coord const& l,
                                                     Box const& box,
                                                     DimsType d);

  static inline Box GetEmptyBox();

  static inline Point GetBoxCenter(Box const& box);

  static Box GetBox(Box const& x, Box const& y);

  template <typename SliceType>
  static Box GetBoxFromSlice(SliceType const V);

  template <typename Range>
  static Box GetBox(Range&& range)
    requires(parlay::is_random_access_range_v<Range>);

  static Box GetBoxFromBoxSeq(BoxSeq const& box_seq);

  template <typename Leaf, typename Interior>
  static Box GetBox(Node* T);

  template <typename Leaf, typename Interior>
  static inline auto RetriveBox(Node const* T)
    requires(HasBox<typename Leaf::AT> && HasBox<typename Interior::AT>);

  //============================================================================
  // SECTION 10: GEOMETRY UTILITIES - CIRCLE OPERATIONS
  // Implementation: base_tree_impl/circle_op.hpp
  //============================================================================

  template <typename CircleType>
  static bool LegalCircle(CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Point const& p, CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Box const& box, CircleType const& cl);

  template <typename CircleType1, typename CircleType2>
  static inline bool CircleWithinCircle(CircleType1 const& a,
                                        CircleType2 const& b);

  template <typename CircleType>
  static inline bool CircleIntersectBox(CircleType const& cl, Box const& box);

  template <typename CircleType>
  static inline bool CircleIntersectCircle(CircleType const& a,
                                           CircleType const& b);

  template <typename CircleType>
  static inline CircleType GetCircle(Box const& box);

  template <typename CircleType>
  static inline CircleType GetCircle(Slice V);

  template <typename CircleType>
  static inline CircleType GetCircle(CircleType const& a, CircleType const& b);

  template <typename CircleType>
  static inline CircleType GetCircle(Point const& p, CircleType const& cl);

  template <typename CircleType>
  static inline CircleType GetCircle(
      parlay::sequence<CircleType> const& circle_seq);

  //============================================================================
  // SECTION 11: GEOMETRY UTILITIES - DISTANCE METRICS
  // Implementation: base_tree_impl/knn_query.hpp
  //============================================================================

  static inline DisType P2PDistanceSquare(Point const& p, Point const& q);

  static inline DisType P2BMinDistanceSquare(Point const& p, Box const& a);

  static inline DisType P2BMaxDistanceSquare(Point const& p, Box const& a);

  static inline double P2CMinDistance(Point const& p, Point const& center,
                                      DisType const r);

  template <typename CircleType>
  static inline double P2CMinDistance(Point const& p, CircleType const& cl);

  static inline DisType InterruptibleDistance(Point const& p, Point const& q,
                                              DisType up);

  //============================================================================
  // SECTION 12: INTERNAL - TREE BUILDING & PARTITIONING
  // Implementation: base_tree_impl/tree_op/build_inner_tree.hpp
  //                 base_tree_impl/points_op.hpp
  //============================================================================

  static inline void SamplePoints(Slice In, Points& arr);

  static inline BucketType FindBucket(Point const& p,
                                      HyperPlaneSeq const& pivots);

  static void Partition(Slice A, Slice B, size_t const n,
                        HyperPlaneSeq const& pivots,
                        parlay::sequence<BallsType>& sums);

  static PointsIter SerialPartition(Slice In, DimsType d);

  template <IsBinaryNode Interior>
  static inline BucketType RetriveTag(Point const& p, NodeTagSeq const& tags);

  template <IsMultiNode Interior>
  static inline BucketType RetriveTag(Point const& p, NodeTagSeq const& tags);

  template <typename Interior>
  static void SeievePoints(Slice A, Slice B, size_t const n,
                           NodeTagSeq const& tags,
                           parlay::sequence<BallsType>& sums,
                           BucketType const tags_num);

  template <typename Leaf, IsBinaryNode Interior, typename Base>
  static Base BuildInnerTree(BucketType idx, HyperPlaneSeq const& pivots,
                             parlay::sequence<Base> const& tree_nodes);

  template <IsMultiNode Interior>
  static Node* BuildInnerTree(BucketType idx, HyperPlaneSeq const& pivots,
                              parlay::sequence<Node*> const& tree_nodes);

  //============================================================================
  // SECTION 13: INTERNAL - NODE OPERATIONS & UPDATES
  // Implementation: base_tree_impl/tree_op/node_op.hpp
  //============================================================================

  template <IsBinaryNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(Node* T, Node* L, Node* R);

  template <IsBinaryNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(Node* T, NodeBox const& L,
                                    NodeBox const& R);

  template <IsMultiNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(
      Node* T, typename Interior::NodeArr const& new_nodes);

  template <typename Leaf, typename Interior>
  static inline typename Interior::ST const& GetSplit(Node const* node)
    requires std::same_as<typename BaseTree::Box, typename Interior::ST>;

  //============================================================================
  // SECTION 14: INTERNAL - TREE TRAVERSAL & EXTRACTION
  // Implementation: base_tree_impl/tree_op/flatten.hpp
  //============================================================================

  template <typename Leaf, IsBinaryNode Interior, typename Range,
            bool granularity = true>
  static void FlattenRec(Node* T, Range Out);

  template <typename Leaf, typename Interior, typename Range,
            bool granularity = true>
  static void FlattenRec(Node* T, Range Out)
    requires(!IsBinaryNode<Interior>);

  template <typename Leaf, typename Range>
  static void ExtractPointsInLeaf(Node* T, Range Out);

  template <typename Leaf, IsMultiNode Interior, typename Range,
            bool granularity = true>
  static void PartialFlatten(Node* T, Range Out, BucketType idx);

  //============================================================================
  // SECTION 15: INTERNAL - TREE TRANSFORMATIONS
  // Implementation: base_tree_impl/tree_op/compress_expand_tree.hpp
  //============================================================================

  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* ExpandMultiNode(typename MN::ST const& split, BucketType idx,
                               BucketType deep,
                               parlay::sequence<Node*> const& tree_nodes);

  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* Expand2Binary(Node* T)
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>;

  template <IsBinaryNode BN, IsMultiNode MN>
  static void CompressBinaryNode(Node* bn_proto,
                                 typename MN::NodeArr& tree_nodes,
                                 typename MN::ST& split, BucketType idx);

  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* Compress2Multi(Node* T)
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>;

  //============================================================================
  // SECTION 16: INTERNAL - IMBALANCE DETECTION & UTILITIES
  // Implementation: base_tree_impl/dimensinality.hpp
  //                 base_tree_impl/box_cut.hpp
  //============================================================================

  static inline size_t GetImbalanceRatio();
  static inline bool ImbalanceNode(size_t const l, size_t const n);
  static inline bool SparcyNode(size_t const l, size_t const n);

  template <SupportsForceParallel Interior, bool granularity>
  inline static bool ForceParallelRecursion(Interior const* T);

 protected:
  //============================================================================
  // SECTION 17: MEMBER VARIABLES
  //============================================================================
  Node* root_ = nullptr;
  parlay::internal::timer timer;
  Box tree_box_;
  size_t delete_node_num_ = 0;
  size_t alloc_node_num_ = 0;
};

}  // namespace psi

//==============================================================================
// IMPLEMENTATION INCLUDES
// Each file implements methods for a specific section above
//==============================================================================

#include "base_tree_impl/box_cut.hpp"                       // Section 16
#include "base_tree_impl/box_op.hpp"                        // Section 9
#include "base_tree_impl/circle_op.hpp"                     // Section 10
#include "base_tree_impl/delete_tree.hpp"                   // Section 4
#include "base_tree_impl/dimensinality.hpp"                 // Section 16
#include "base_tree_impl/inner_tree/inner_tree.hpp"         // Section 3, 12
#include "base_tree_impl/knn_query.hpp"                     // Section 5.1, 11
#include "base_tree_impl/points_op.hpp"                     // Section 12
#include "base_tree_impl/range_query.hpp"                   // Section 5.2
#include "base_tree_impl/tree_op/build_inner_tree.hpp"      // Section 4, 12
#include "base_tree_impl/tree_op/compress_expand_tree.hpp"  // Section 15
#include "base_tree_impl/tree_op/flatten.hpp"               // Section 14
#include "base_tree_impl/tree_op/leaf_op.hpp"               // Section 6
#include "base_tree_impl/tree_op/node_op.hpp"               // Section 13
#include "base_tree_impl/tree_op/rebuild.hpp"               // Section 6
#include "base_tree_impl/validation.hpp"                    // Section 7, 8

#endif  // PSI_BASE_TREE_H_
