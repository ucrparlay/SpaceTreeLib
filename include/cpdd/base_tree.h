#pragma once

#include <sys/types.h>

#include <cstdint>
#include <type_traits>

#include "dependence/comparator.h"
#include "dependence/concepts.h"
#include "dependence/loggers.h"
#include "dependence/search_container.h"
#include "dependence/tree_node.h"

namespace cpdd {

#define LOG std::cout
#define ENDL std::endl << std::flush

template <typename Point, typename DerivedTree, uint_fast8_t kBDO = 6>
class BaseTree {
 public:
  // NOTE: when kBDO >= 8, the # bucket is 255, total skeleton nodes >= 255*2
  using BucketType =
      std::conditional_t<(kBDO > 7), uint_fast16_t, uint_fast8_t>;
  using BallsType = uint_fast32_t;
  using DimsType = uint_fast8_t;
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
  using Circle = std::pair<Point, Coord>;

  using NodeBox = std::pair<Node*, Box>;
  using NodeBoxSeq = parlay::sequence<NodeBox>;
  using NodeTag = std::pair<Node*, uint_fast8_t>;
  using NodeTagSeq = parlay::sequence<NodeTag>;

  // NOTE: Const variables
  // NOTE: uint32t handle up to 4e9 at least
  // WARN: bucket num should smaller than 1<<8 to handle type overflow
  static constexpr DimsType const kDim = std::tuple_size_v<Coords>;
  static constexpr BucketType const kBuildDepthOnce = kBDO;
  static constexpr BucketType const kPivotNum = (1 << kBuildDepthOnce) - 1;
  static constexpr BucketType const kBucketNum = 1 << kBuildDepthOnce;

  // NOTE: tree structure
  static constexpr uint_fast8_t const kLeaveWrap = 32;
  static constexpr uint_fast8_t const kThinLeaveWrap = 24;
  static constexpr uint_fast16_t const kSerialBuildCutoff = 1 << 10;

  // NOTE: block param in Partition
  static constexpr uint_fast8_t const kLog2Base = 10;
  static constexpr uint_fast16_t const kBlockSize = 1 << kLog2Base;

  // NOTE: reconstruct weight threshold
  static constexpr uint_fast8_t const kInbalanceRatio = 30;

  // NOTE: get the imbalance ratio
  static inline size_t GetImbalanceRatio();
  static inline bool ImbalanceNode(size_t const l, size_t const n);
  static inline bool SparcyNode(size_t const l, size_t const n);

  // NOTE: array based inner tree for batch insertion and deletion
  template <typename Leaf, typename Interior>
  struct InnerTree;

  struct BoxCut;

  // NOTE: Box operations
  static inline bool LegalBox(Box const& bx);
  static inline bool WithinBox(Box const& a, Box const& b);
  static inline bool WithinBox(Point const& p, Box const& bx);
  static inline bool BoxIntersectBox(Box const& a, Box const& b);
  static inline Box GetEmptyBox();
  static Box GetBox(Box const& x, Box const& y);
  static Box GetBox(Slice V);
  template <typename Leaf, typename Interior>
  static Box GetBox(Node* T);
  static Box GetBox(BoxSeq const& box_seq);

  static inline bool WithinCircle(Box const& bx, Circle const& cl);
  static inline bool WithinCircle(Point const& p, Circle const& cl);
  static inline bool CircleIntersectBox(Circle const& cl, Box const& bx);

  // NOTE: build tree
  static inline void SamplePoints(Slice In, Points& arr);

  static inline BucketType FindBucket(Point const& p,
                                      HyperPlaneSeq const& pivots);

  template <IsBinaryNode Interior>
  static inline void UpdateInterior(Node* T, Node* L, Node* R);

  template <IsBinaryNode Interior>
  static inline void UpdateInterior(Node* T, NodeBox const& L,
                                    NodeBox const& R);

  template <IsMultiNode Interior>
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

  template <typename Leaf, typename Interior, bool granularity = true,
            typename... Args>
  Node* RebuildSingleTree(Node* T, Args&&... args);

  // virtual Node* BuildRecursiveWrapper(Slice In, Slice Out, const Box& bx,
  //                                     DimsType dim) = 0;

  template <IsBinaryNode Interior>
  static Node* BuildInnerTree(BucketType idx, HyperPlaneSeq& pivots,
                              parlay::sequence<Node*>& tree_nodes);

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

  // NOTE: KNN query stuffs
  static inline Coord P2PDistance(Point const& p, Point const& q);

  static inline Coord P2BMinDistance(Point const& p, Box const& a);

  static inline Coord P2BMaxDistance(Point const& p, Box const& a);

  static inline Coord InterruptibleDistance(Point const& p, Point const& q,
                                            Coord up);

  // NOTE: searech knn in the leaf
  template <typename Leaf, typename Range>
  static void KNNLeaf(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                      Box const& bx);

  // NOTE: search knn in the binary node
  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void KNNBinary(Node* T, Point const& q,
                        kBoundedQueue<Point, Range>& bq, Box const& bx,
                        KNNLogger& logger);

  // NOTE: search knn in the expanded multi node
  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMultiExpand(Node* T, Point const& q, DimsType dim,
                             BucketType idx, kBoundedQueue<Point, Range>& bq,
                             Box const& bx, KNNLogger& logger);

  // NOTE: search knn in the multi node
  template <typename Leaf, IsMultiNode Interior, typename Range>
  static void KNNMulti(Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
                       Box const& bx, KNNLogger& logger);

  // NOTE: range count stuffs
  template <typename Leaf>
  static size_t RangeCountRectangleLeaf(Node* T, Box const& query_box,
                                        Box const& node_box);

  template <typename Leaf, IsBinaryNode Interior>
  static size_t RangeCountRectangle(Node* T, Box const& query_box,
                                    Box const& node_box,
                                    RangeQueryLogger& logger);

  template <typename Leaf, IsMultiNode Interior>
  static size_t RangeCountRectangle(Node* T, Box const& query_box,
                                    Box const& node_box, DimsType dim,
                                    BucketType idx, RangeQueryLogger& logger);

  template <typename Leaf, IsBinaryNode Interior>
  static size_t RangeCountRadius(Node* T, Circle const& cl,
                                 Box const& node_box);

  // NOTE: range query stuffs
  template <typename Leaf, typename Range>
  static void RangeQueryLeaf(Node* T, Range Out, size_t& s,
                             Box const& query_box, Box const& node_box);

  template <typename Leaf, IsBinaryNode Interior, typename Range>
  static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                        Box const& query_box,
                                        Box const& node_box,
                                        RangeQueryLogger& logger);

  template <typename Leaf, IsMultiNode Interior>
  static size_t RangeCountRadius(Node* T, Circle const& cl,
                                 Box const& node_box);

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

  template <typename Leaf, IsMultiNode Interior, typename Range,
            bool granularity = true>
  static void FlattenRec(Node* T, Range Out);

  template <typename Leaf, IsMultiNode Interior, typename Range,
            bool granularity = true>
  static void PartialFlatten(Node* T, Range Out, BucketType idx);

  template <IsBinaryNode BN, IsMultiNode MN>
  static Node* ExpandMultiNode(typename MN::ST const& split, BucketType idx,
                               BucketType deep,
                               parlay::sequence<Node*> const& tree_nodes);

  template <IsBinaryNode BN, IsMultiNode MN>
    requires std::same_as<typename BN::ST, typename MN::ST::value_type>
  static Node* Expand2Binary(Node* T);

  // NOTE: validations
  template <typename Leaf, typename Interior>
  Box CheckBox(Node* T, Box const& box);

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

  // NOTE: interfaces
  inline void SetRoot(Node* root) { this->root_ = root; }

  inline Node* GetRoot() { return this->root_; }

  inline Box GetRootBox() { return this->tree_box_; }

 protected:
  Node* root_ = nullptr;
  parlay::internal::timer timer;
  Box tree_box_;
  size_t delete_node_num_ = 0;
  size_t alloc_node_num_ = 0;
};

}  // namespace cpdd

#include "base_tree_impl/box_cut.hpp"
#include "base_tree_impl/box_op.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/inner_tree.hpp"
#include "base_tree_impl/knn_query.hpp"
#include "base_tree_impl/points_op.hpp"
#include "base_tree_impl/range_query.hpp"
#include "base_tree_impl/tree_op.hpp"
#include "base_tree_impl/validation.hpp"
