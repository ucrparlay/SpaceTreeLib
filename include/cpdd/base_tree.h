#pragma once

#include <sys/types.h>
#include <cstdint>
#include "dependence/comparator.h"
#include "dependence/tree_node.h"
#include "dependence/search_container.h"

namespace cpdd {

#define LOG  std::cout
#define ENDL std::endl << std::flush

template<typename Point, uint8_t kBDO = 6>
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
    /*using Splitter = std::pair<Coord, DimsType>;*/
    /*using SplitterSeq = parlay::sequence<Splitter>;*/
    using Box = std::pair<Point, Point>;
    using BoxSeq = parlay::sequence<Box>;
    using Circle = std::pair<Point, Coord>;

    using NodeBox = std::pair<Node*, Box>;
    using NodeTag = std::pair<Node*, uint_fast8_t>;
    using NodeTagSeq = parlay::sequence<NodeTag>;
    using TagNodes = parlay::sequence<BallsType>;

    //@ Const variables
    //@ uint32t handle up to 4e9 at least
    //! bucket num should smaller than 1<<8 to handle type overflow

    // TODO: wrap the variables using std::getenv()
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
    inline size_t GetImbalanceRatio();
    inline bool ImbalanceNode(const size_t l, const size_t n);

    // NOTE: array based inner tree for batch insertion and deletion
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

    // NOTE: dimensionality
    template<typename SplitRule>
    inline DimsType PickRebuildDim(const Node* T, const DimsType d,
                                   const DimsType DIM);

    // NOTE: build tree
    static inline void SamplePoints(Slice In, Points& arr);

    template<typename SplitterSeq>
    static inline BucketType FindBucket(const Point& p,
                                        const SplitterSeq& pivots);

    template<typename SplitterSeq>
    static void Partition(Slice A, Slice B, const size_t n,
                          const SplitterSeq& pivots,
                          parlay::sequence<BallsType>& sums);

    template<typename Interior, typename SplitterSeq>
    static Node* BuildInnerTree(BucketType idx, SplitterSeq& pivots,
                                parlay::sequence<Node*>& tree_nodes);

    PointsIter SerialPartition(Slice In, DimsType d);

    // NOTE: delete tree
    virtual void DeleteTree() = 0;

    template<typename leaf_node_type, typename interior_node_type>
    void DeleteTreeWrapper();

    template<typename Leaf, typename Interior>
    static void DeleteTreeRecursive(Node* T, bool granularity = true);

    // NOTE: KNN query stuffs
    static inline Coord P2PDistance(const Point& p, const Point& q,
                                    const DimsType DIM);

    static inline Coord P2BMinDistance(const Point& p, const Box& a,
                                       const DimsType DIM);

    static inline Coord P2BMaxDistance(const Point& p, const Box& a,
                                       const DimsType DIM);

    static inline Coord InterruptibleDistance(const Point& p, const Point& q,
                                              Coord up, DimsType DIM);

    template<typename Leaf, typename Interior, typename StoreType>
    static void KNNRec(Node* T, const Point& q, const DimsType DIM,
                       kBoundedQueue<Point, StoreType>& bq, const Box& bx,
                       size_t& vis_node_num);

    // NOTE: range count stuffs
    template<typename Leaf, typename Interior>
    static size_t RangeCountRectangle(Node* T, const Box& query_box,
                                      const Box& node_box, size_t& vis_leaf_num,
                                      size_t& vis_inter_num);

    template<typename Leaf, typename Interior>
    static size_t RangeCountRadius(Node* T, const Circle& cl,
                                   const Box& node_box);

    // NOTE: range query stuffs
    template<typename Leaf, typename Interior, typename StoreType>
    static void RangeQuerySerialRecursive(Node* T, StoreType Out, size_t& s,
                                          const Box& query_box,
                                          const Box& node_box);

    // NOTE: utility
    template<typename Leaf, typename Interior, typename Range>
    static void FlattenRec(Node* T, Range Out, bool granularity = true);

    // NOTE: validations
    template<typename Leaf, typename Interior>
    bool CheckBox(Node* T, const Box& bx);

    template<typename Leaf, typename Interior>
    size_t CheckSize(Node* T);

    template<typename Leaf, typename Interior>
    void CheckTreeSameSequential(Node* T, int dim, const int& DIM);

    template<typename Leaf, typename Interior>
    void Validate(const DimsType DIM);

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

    //@ kdtree interfaces
    inline void SetRoot(Node* _root) { this->root_ = _root; }

    inline Node* GetRoot() { return this->root_; }

    inline Box GetRootBox() { return this->tree_box_; }

 protected:
    Node* root_ = nullptr;
    parlay::internal::timer timer;
    Box tree_box_;
};

}  // namespace cpdd

#include "base_tree_impl/box_op.hpp"
#include "base_tree_impl/validation.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/points_op.hpp"
#include "base_tree_impl/tree_op.hpp"
#include "base_tree_impl/query_op.hpp"
