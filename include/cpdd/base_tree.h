#pragma once

#include <sys/types.h>
#include <cstdint>
#include "dependence/comparator.h"
#include "dependence/tree_node.h"
#include "dependence/search_container.h"
#include "dependence/concepts.h"
#include "dependence/loggers.h"

namespace cpdd {

#define LOG  std::cout
#define ENDL std::endl << std::flush

template<typename Point, uint_fast8_t kBDO = 6>
class BaseTree {
 public:
    using BucketType = uint_fast8_t;
    using BallsType = uint_fast32_t;
    using DimsType = uint_fast8_t;

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
    using NodeTag = std::pair<Node*, uint_fast8_t>;
    using NodeTagSeq = parlay::sequence<NodeTag>;
    using TagNodes = parlay::sequence<BallsType>;

    // NOTE: Const variables
    // NOTE: uint32t handle up to 4e9 at least
    // WARN: bucket num should smaller than 1<<8 to handle type overflow
    static const constexpr DimsType kDim = std::tuple_size_v<Coords>;
    static const constexpr BucketType kBuildDepthOnce = kBDO;
    static const constexpr BucketType kPivotNum = (1 << kBuildDepthOnce) - 1;
    static const constexpr BucketType kBucketNum = 1 << kBuildDepthOnce;

    // NOTE: tree structure
    static const constexpr uint_fast8_t kLeaveWrap = 32;
    static const constexpr uint_fast8_t kThinLeaveWrap = 24;
    static const constexpr uint_fast16_t kSerialBuildCutoff = 1 << 10;

    // NOTE: block param in Partition
    static const constexpr uint_fast8_t kLog2Base = 10;
    static const constexpr uint_fast16_t kBlockSize = 1 << kLog2Base;

    // NOTE: reconstruct weight threshold
    static const constexpr uint_fast8_t kInbalanceRatio = 30;

    // NOTE: get the imbalance ratio
    static inline size_t GetImbalanceRatio();
    static inline bool ImbalanceNode(const size_t l, const size_t n);

    // NOTE: array based inner tree for batch insertion and deletion
    template<typename Leaf, typename Interior>
    struct InnerTree;

    // NOTE: Box operations
    static inline bool LegalBox(const Box& bx);
    static inline bool WithinBox(const Box& a, const Box& b);
    static inline bool WithinBox(const Point& p, const Box& bx);
    static inline bool BoxIntersectBox(const Box& a, const Box& b);
    static inline Box GetEmptyBox();
    static Box GetBox(const Box& x, const Box& y);
    static Box GetBox(Slice V);
    static Box GetBox(Node* T);

    static inline bool WithinCircle(const Box& bx, const Circle& cl);
    static inline bool WithinCircle(const Point& p, const Circle& cl);
    static inline bool CircleIntersectBox(const Circle& cl, const Box& bx);

    // NOTE: build tree
    static inline void SamplePoints(Slice In, Points& arr);

    static inline BucketType FindBucket(const Point& p,
                                        const HyperPlaneSeq& pivots);

    template<IsBinaryNode Interior>
    static inline void UpdateInterior(Node* T, Node* L, Node* R);

    template<IsMultiNode Interior>
    static inline void UpdateInterior(Node* T,
                                      typename Interior::Nodes& new_nodes);

    template<typename Leaf, typename Interior, bool granularity = true>
    static void PrepareRebuild(Node* T, Slice In, Points& wx, Points& wo);

    static void Partition(Slice A, Slice B, const size_t n,
                          const HyperPlaneSeq& pivots,
                          parlay::sequence<BallsType>& sums);

    // NOTE: batch insert
    template<typename Leaf>
    static Node* InsertPoints2Leaf(Node* T, Slice In);

    template<IsBinaryNode Interior>
    static inline BucketType RetriveTag(const Point& p, const NodeTagSeq& tags);

    template<IsMultiNode Interior>
    static inline BucketType RetriveTag(const Point& p, const NodeTagSeq& tags);

    template<typename Interior>
    static void SeievePoints(Slice A, Slice B, const size_t n,
                             const NodeTagSeq& tags,
                             parlay::sequence<BallsType>& sums,
                             const BucketType tags_num);

    template<IsBinaryNode Interior>
    static Node* BuildInnerTree(BucketType idx, HyperPlaneSeq& pivots,
                                parlay::sequence<Node*>& tree_nodes);

    static PointsIter SerialPartition(Slice In, DimsType d);

    // NOTE: delete tree
    template<SupportsForceParallel Interior, bool granularity>
    inline static bool ForceParallelRecursion(Interior* T);

    virtual void DeleteTree() = 0;

    template<typename Leaf, typename Interior>
    void DeleteTreeWrapper();

    template<typename Leaf, IsBinaryNode Interior, bool granularity = true>
    static void DeleteTreeRecursive(Node* T);

    template<typename Leaf, IsMultiNode Interior, bool granularity = true>
    static void DeleteTreeRecursive(Node* T);

    // NOTE: KNN query stuffs
    static inline Coord P2PDistance(const Point& p, const Point& q);

    static inline Coord P2BMinDistance(const Point& p, const Box& a);

    static inline Coord P2BMaxDistance(const Point& p, const Box& a);

    static inline Coord InterruptibleDistance(const Point& p, const Point& q,
                                              Coord up);

    // NOTE: searech knn in the leaf
    template<typename Leaf, typename Range>
    static void KNNLeaf(Node* T, const Point& q,
                        kBoundedQueue<Point, Range>& bq, const Box& bx);

    // NOTE: search knn in the binary node
    template<typename Leaf, IsBinaryNode Interior, typename Range>
    static void KNNBinary(Node* T, const Point& q,
                          kBoundedQueue<Point, Range>& bq, const Box& bx,
                          KNNLogger& logger);

    // NOTE: search knn in the expanded multi node
    template<typename Leaf, IsMultiNode Interior, typename Range>
    static void KNNMultiExpand(Node* T, const Point& q, DimsType dim,
                               BucketType idx, kBoundedQueue<Point, Range>& bq,
                               const Box& bx, KNNLogger& logger);

    // NOTE: search knn in the multi node
    template<typename Leaf, IsMultiNode Interior, typename Range>
    static void KNNMulti(Node* T, const Point& q,
                         kBoundedQueue<Point, Range>& bq, const Box& bx,
                         KNNLogger& logger);

    // NOTE: range count stuffs
    template<typename Leaf>
    static size_t RangeCountRectangleLeaf(Node* T, const Box& query_box,
                                          const Box& node_box);

    template<typename Leaf, IsBinaryNode Interior>
    static size_t RangeCountRectangle(Node* T, const Box& query_box,
                                      const Box& node_box,
                                      RangeQueryLogger& logger);

    template<typename Leaf, IsMultiNode Interior>
    static size_t RangeCountRectangle(Node* T, const Box& query_box,
                                      const Box& node_box, DimsType dim,
                                      BucketType idx, RangeQueryLogger& logger);

    template<typename Leaf, IsBinaryNode Interior>
    static size_t RangeCountRadius(Node* T, const Circle& cl,
                                   const Box& node_box);

    // NOTE: range query stuffs
    template<typename Leaf, typename Range>
    static void RangeQueryLeaf(Node* T, Range Out, size_t& s,
                               const Box& query_box, const Box& node_box);

    template<typename Leaf, IsBinaryNode Interior, typename Range>
    static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                          const Box& query_box,
                                          const Box& node_box,
                                          RangeQueryLogger& logger);

    template<typename Leaf, IsMultiNode Interior>
    static size_t RangeCountRadius(Node* T, const Circle& cl,
                                   const Box& node_box);

    // NOTE: range query stuffs
    template<typename Leaf, IsMultiNode Interior, typename Range>
    static void RangeQuerySerialRecursive(Node* T, Range Out, size_t& s,
                                          const Box& query_box,
                                          const Box& node_box, DimsType dim,
                                          BucketType idx,
                                          RangeQueryLogger& logger);

    // NOTE: utility
    // TODO: better evaluate the parallel recursion function
    template<typename Leaf, IsBinaryNode Interior, typename Range,
             bool granularity = true>
    static void FlattenRec(Node* T, Range Out);

    template<typename Leaf, IsMultiNode Interior, typename Range,
             bool granularity = true>
    static void FlattenRec(Node* T, Range Out);

    template<typename Leaf, IsMultiNode Interior, typename Range,
             bool granularity = true>
    static void PartialFlatten(Node* T, Range Out, BucketType idx);

    template<IsBinaryNode BN, IsMultiNode MN>
    static Node* ExpandMultiNode(const typename MN::ST& split, BucketType idx,
                                 BucketType deep,
                                 const parlay::sequence<Node*>& tree_nodes);

    template<IsBinaryNode BN, IsMultiNode MN>
        requires std::same_as<typename BN::ST, typename MN::ST::value_type>
    static Node* Expand2Binary(Node* T);

    // NOTE: validations
    template<typename Leaf, typename Interior>
    bool CheckBox(Node* T, const Box& bx);

    template<typename Leaf, typename Interior>
    size_t CheckSize(Node* T);

    template<typename Leaf, typename Interior>
    void CheckTreeSameSequential(Node* T, int dim);

    template<typename Leaf, typename Interior, typename SplitRule>
    void Validate();

    template<typename Leaf, typename Interior>
    size_t GetTreeHeight();

    template<typename Leaf, typename Interior>
    size_t GetMaxTreeDepth(Node* T, size_t deep);

    template<typename Leaf, typename Interior>
    double GetAveTreeHeight();

    template<typename Leaf, typename Interior>
    size_t CountTreeNodesNum(Node* T);

    template<typename Leaf, typename Interior>
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

#include "base_tree_impl/box_op.hpp"
#include "base_tree_impl/validation.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/inner_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/points_op.hpp"
#include "base_tree_impl/tree_op.hpp"
#include "base_tree_impl/knn_query.hpp"
#include "base_tree_impl/range_query.hpp"
