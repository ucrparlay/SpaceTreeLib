#ifndef PSPT_BASE_TREE_H_
#define PSPT_BASE_TREE_H_

#include <sys/types.h>

#include <cstdint>
#include <type_traits>

#include "dependence/comparator.h"
#include "dependence/concepts.h"
#include "dependence/loggers.h"
#include "dependence/search_container.h"
#include "dependence/tree_node.h"

namespace pspt {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class BaseTree {
 public:
  // NOTE: when kSkHeight >= 8, the # bucket is 255, total skeleton nodes >=
  // 255*2
  using BucketType =
      std::conditional_t<(kSkHeight > 7), uint_fast16_t, uint_fast8_t>;
  using BallsType = uint_fast32_t;
  using DimsType = uint_fast8_t;
  // using DepthType = uint_fast8_t;
  using DepthType = int;
  using BucketSeq = parlay::sequence<BucketType>;
  using BallSeq = parlay::sequence<BallsType>;

  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Num = Num_Comparator<Coord>;
  using Slice = parlay::slice<Point*, Point*>;
  using Points = parlay::sequence<Point>;
  using PointsIter = typename parlay::sequence<Point>::iterator;
  using HyperPlane = std::pair<Coord, DimsType>;
  using HyperPlaneSeq = parlay::sequence<HyperPlane>;
  using Box = std::pair<Point, Point>;
  using BoxSeq = parlay::sequence<Box>;

  using NodeBoolean = std::pair<Node*, bool>;
  using NodeBox = std::pair<Node*, Box>;
  using NodeBoxSeq = parlay::sequence<NodeBox>;
  using NodeTag = std::pair<Node*, uint_fast8_t>;
  using NodeTagSeq = parlay::sequence<NodeTag>;

  // NOTE: Const variables
  // NOTE: uint32t handle up to 4e9 at least
  // WARN: bucket num should smaller than 1<<8 to handle type overflow
  static constexpr DimsType const kDim = std::tuple_size_v<Coords>;
  static constexpr BucketType const kBuildDepthOnce = kSkHeight;
  static constexpr BucketType const kPivotNum = (1 << kBuildDepthOnce) - 1;
  static constexpr BucketType const kBucketNum = 1 << kBuildDepthOnce;

  // NOTE: tree structure
  // static constexpr uint_fast8_t const kLeaveWrap = 32;
  // static constexpr uint_fast8_t const kThinLeaveWrap = 24;
  // static constexpr uint_fast8_t const kSlimLeaveWrap = 8;
  static constexpr uint_fast8_t const kLeaveWrap = 8;
  static constexpr uint_fast8_t const kThinLeaveWrap = 4;
  static constexpr uint_fast8_t const kSlimLeaveWrap = 2;
  static constexpr uint_fast16_t const kSerialBuildCutoff = 1 << 10;

  // NOTE: block param in Partition
  static constexpr uint_fast8_t const kLog2Base = 10;
  static constexpr uint_fast16_t const kBlockSize = 1 << kLog2Base;

  // NOTE: reconstruct weight threshold
  static constexpr uint_fast8_t const kInbalanceRatio = kImbaRatio;

  // NOTE: array based inner tree for batch insertion and deletion
  template <typename Leaf, typename Interior>
  struct InnerTree;

  // NOTE: compute the bounding box on the fly
  struct BoxCut;

  // NOTE: get the imbalance ratio
  static inline size_t GetImbalanceRatio();
  static inline bool ImbalanceNode(size_t const l, size_t const n);
  static inline bool SparcyNode(size_t const l, size_t const n);

  // NOTE: Box operations
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
  static Box GetBox(Slice V);
  static Box GetBox(BoxSeq const& box_seq);
  template <typename Leaf, typename Interior>
  static Box GetBox(Node* T);

  // NOTE: Circle operations
  struct NormalCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return radius; }

    Coord GetRadiusSquare() const { return radius * radius; }

    static Coord ComputeRadius(double const& r) {
      return static_cast<Coord>(r);
    }

    Point center;
    Coord radius;
  };

  struct CoverCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return static_cast<Coord>(1 << level); }

    Coord GetRadiusSquare() const { return GetRadius() * GetRadius(); }

    static DepthType ComputeRadius(double const& r) {
      return static_cast<DepthType>(std::ceil(std::log2(r)));
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

  template <typename CircleType>
  static bool LegalCircle(CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Point const& p, CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Box const& box, CircleType const& cl);

  template <typename CircleType>
  static inline bool CircleWithinCircle(CircleType const& a,
                                        CircleType const& b);

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

  // NOTE: build tree
  static inline void SamplePoints(Slice In, Points& arr);

  static inline BucketType FindBucket(Point const& p,
                                      HyperPlaneSeq const& pivots);

  template <IsBinaryNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(Node* T, Node* L, Node* R);

  template <IsBinaryNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(Node* T, NodeBox const& L,
                                    NodeBox const& R);

  template <IsMultiNode Interior, bool UpdateParFlag = true>
  static inline void UpdateInterior(
      Node* T, typename Interior::NodeArr const& new_nodes);

  template <typename Leaf, typename Interior, bool granularity = true>
  static void PrepareRebuild(Node* T, Slice In, Points& wx, Points& wo);

  template <typename Leaf, typename Interior, bool granularity = true>
  static void PrepareRebuild(Node* T, Points& wx, Points& wo);

  static void Partition(Slice A, Slice B, size_t const n,
                        HyperPlaneSeq const& pivots,
                        parlay::sequence<BallsType>& sums);

  // NOTE: batch insert
  template <typename Leaf>
  static Node* InsertPoints2Leaf(Node* T, Slice In);

  template <typename Leaf, typename RT>
  static RT DeletePoints4Leaf(Node* T, Slice In);

  template <typename Leaf, typename RT>
  static RT DiffPoints4Leaf(Node* T, Slice In);

  template <IsBinaryNode Interior>
  static inline BucketType RetriveTag(Point const& p, NodeTagSeq const& tags);

  template <IsMultiNode Interior>
  static inline BucketType RetriveTag(Point const& p, NodeTagSeq const& tags);

  template <typename Interior>
  static void SeievePoints(Slice A, Slice B, size_t const n,
                           NodeTagSeq const& tags,
                           parlay::sequence<BallsType>& sums,
                           BucketType const tags_num);

  template <typename Leaf, typename Interior, typename... Args>
  Node* RebuildWithInsert(Node* T, Slice In, Args&&... args);

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

  template <typename Leaf, IsBinaryNode Interior, typename Base>
  static Base BuildInnerTree(BucketType idx, HyperPlaneSeq const& pivots,
                             parlay::sequence<Base> const& tree_nodes);

  template <IsMultiNode Interior>
  static Node* BuildInnerTree(BucketType idx, HyperPlaneSeq const& pivots,
                              parlay::sequence<Node*> const& tree_nodes);

  static PointsIter SerialPartition(Slice In, DimsType d);

  // NOTE: delete tree
  template <SupportsForceParallel Interior, bool granularity>
  inline static bool ForceParallelRecursion(Interior const* T);

  constexpr virtual void DeleteTree() = 0;

  template <typename Leaf, typename Interior>
  void DeleteTreeWrapper();

  template <typename Leaf, IsBinaryNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  template <typename Leaf, IsMultiNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  template <typename Leaf, IsDynamicNode Interior, bool granularity = true>
  static void DeleteTreeRecursive(Node* T);

  // NOTE: KNN query stuffs
  static inline Coord P2PDistanceSquare(Point const& p, Point const& q);

  static inline Coord P2BMinDistance(Point const& p, Box const& a);

  static inline Coord P2BMaxDistance(Point const& p, Box const& a);

  static inline Coord P2CMinDistance(Point const& p, Point const& center,
                                     Coord const r);

  template <typename CircleType>
  static inline Coord P2CMinDistance(Point const& p, CircleType const& cl);

  static inline Coord InterruptibleDistance(Point const& p, Point const& q,
                                            Coord up);

  // NOTE: get the split for a node
  template <typename Leaf, typename Interior>
  static inline typename Interior::ST const& GetSplit(Node const* node)
    requires std::same_as<typename BaseTree::Box, typename Interior::ST>;

  // NOTE: searech knn in the leaf
  template <typename Leaf, typename Range>
  static void KNNLeaf(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq);

  // NOTE: search knn in the binary node
  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void KNNBinary(Node* T, Point const& q,
                        kBoundedQueue<Point, Range>& bq, Box const& node_box,
                        KNNLogger& logger)
    requires std::same_as<typename Interior::ST, HyperPlane>;

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void KNNBinary(Node* T, Point const& q,
                        kBoundedQueue<Point, Range>& bq, KNNLogger& logger)
    requires std::same_as<typename Interior::ST, Box>;

  // NOTE: search knn in the expanded multi node
  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMultiExpand(Node* T, Point const& q, DimsType dim,
                             BucketType idx, kBoundedQueue<Point, Range>& bq,
                             Box const& node_box, KNNLogger& logger);

  // NOTE: search knn in the multi node
  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMulti(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                       Box const& node_box, KNNLogger& logger);

  // NOTE: search knn in the mix of binary and multi node
  template <typename Leaf, IsBinaryNode BN, IsMultiNode MN, typename Range>
  static void KNNMix(Node* T, Point const& q, DimsType dim, BucketType idx,
                     kBoundedQueue<Point, Range>& bq, Box const& node_box,
                     KNNLogger& logger);

  // NOTE: search knn in the cover node
  // TODO: change type of interior
  template <typename Leaf, IsDynamicNode Interior, typename Range>
  static void KNNCover(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                       KNNLogger& logger);

  // NOTE: range count stuffs
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

  // template <typename Leaf, IsBinaryNode Interior>
  // static size_t RangeCountRadius(Node* T, NormalCircle const& cl,
  //                                Box const& node_box);

  // NOTE: range query stuffs
  template <typename Leaf, typename Range>
  static void RangeQueryLeaf(Node* T, Range Out, size_t& s,
                             Box const& query_box);

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                        Box const& query_box,
                                        Box const& node_box,
                                        RangeQueryLogger& logger);

  // template <typename Leaf, IsMultiNode Interior>
  // static size_t RangeCountRadius(Node* T, NormalCircle const& cl,
  //                                Box const& node_box);

  // NOTE: range query stuffs
  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                        Box const& query_box,
                                        Box const& node_box, DimsType dim,
                                        BucketType idx,
                                        RangeQueryLogger& logger);

  // NOTE: utility
  // TODO: better evaluate the parallel recursion function
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

  // NOTE: expand a multi node into a binary node
  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* ExpandMultiNode(typename MN::ST const& split, BucketType idx,
                               BucketType deep,
                               parlay::sequence<Node*> const& tree_nodes);

  // NOTE: expand a multi-way tree into a binary tree
  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* Expand2Binary(Node* T)
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>;

  // NOTE: compress several level of bianry node into a multi node
  template <IsBinaryNode BN, IsMultiNode MN>
  static void CompressBinaryNode(Node* bn_proto,
                                 typename MN::NodeArr& tree_nodes,
                                 typename MN::ST& split, BucketType idx);

  // NOTE: compress a binary tree into a multi-way tree
  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* Compress2Multi(Node* T)
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>;

  // NOTE: validations
  template <typename Leaf, typename Interior>
  Box CheckBox(Node* T, Box const& box);

  template <typename Leaf, typename Interior>
  typename Interior::CircleType CheckCover(
      Node* T, typename Interior::CircleType const& level_cover_circle);

  template <typename Leaf, typename Interior>
  static size_t CheckSize(Node* T);

  template <typename Leaf, typename Interior>
  void CheckTreeSameSequential(Node* T, int dim);

  template <typename Leaf, typename Interior, typename SplitRule>
  void Validate();

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

  // NOTE: param interfaces
  void SetRoot(Node* root) { this->root_ = root; }

  Node* GetRoot() { return this->root_; }

  Box GetRootBox() { return this->tree_box_; }

  consteval static auto GetBuildDepthOnce() {
    return static_cast<int>(kBuildDepthOnce);
  }

 protected:
  Node* root_ = nullptr;
  parlay::internal::timer timer;
  Box tree_box_;
  size_t delete_node_num_ = 0;
  size_t alloc_node_num_ = 0;
};

}  // namespace pspt

#include "base_tree_impl/box_cut.hpp"
#include "base_tree_impl/box_op.hpp"
#include "base_tree_impl/circle_op.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/inner_tree.hpp"
#include "base_tree_impl/knn_query.hpp"
#include "base_tree_impl/points_op.hpp"
#include "base_tree_impl/range_query.hpp"
#include "base_tree_impl/tree_op/build_inner_tree.hpp"
#include "base_tree_impl/tree_op/compress_expand_tree.hpp"
#include "base_tree_impl/tree_op/flatten.hpp"
#include "base_tree_impl/tree_op/leaf_op.hpp"
#include "base_tree_impl/tree_op/node_op.hpp"
#include "base_tree_impl/tree_op/rebuild.hpp"
#include "base_tree_impl/validation.hpp"

#endif  // PSPT_BASE_TREE_H_
