#pragma once

#include "dependence/comparator.h"
#include "dependence/tree_node.h"
#include "dependence/search_container.h"

namespace cpdd {

#define LOG  std::cout
#define ENDL std::endl << std::flush

template<typename Point>
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
    static constexpr BucketType kBuildDepthOnce = 6;  //* last layer is leaf
    static constexpr BucketType kPivotNum = (1 << kBuildDepthOnce) - 1;
    static constexpr BucketType kBucketNum = 1 << kBuildDepthOnce;
    // NOTE: tree structure
    static constexpr uint_fast8_t kLeaveWrap = 32;
    static constexpr uint_fast8_t kThinLeaveWrap = 24;
    static constexpr uint_fast16_t kSerialBuildCutoff = 1 << 10;
    // NOTE: block param in Partition
    static constexpr uint_fast8_t kLog2Base = 10;
    static constexpr uint_fast16_t kBlockSize = 1 << kLog2Base;
    // NOTE: reconstruct weight threshold
    static constexpr uint_fast8_t kInbalanceRatio = 30;

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

    template<typename InterType, typename SplitterSeq>
    static Node* BuildInnerTree(BucketType idx, SplitterSeq& pivots,
                                parlay::sequence<Node*>& tree_nodes);

    PointsIter SerialPartition(Slice In, DimsType d);

    // NOTE: delete tree
    virtual void DeleteTree() = 0;

    template<typename leaf_node_type, typename interior_node_type>
    void DeleteTreeWrapper();

    template<typename leaf_node_type, typename interior_node_type>
    static void DeleteTreeRecursive(Node* T, bool granularity = true);

    // NOTE: query stuffs
    static inline Coord P2PDistance(const Point& p, const Point& q,
                                    const DimsType DIM);

    static inline Coord P2BMinDistance(const Point& p, const Box& a,
                                       const DimsType DIM);

    static inline Coord P2BMaxDistance(const Point& p, const Box& a,
                                       const DimsType DIM);

    static inline Coord InterruptibleDistance(const Point& p, const Point& q,
                                              Coord up, DimsType DIM);

    template<typename LeafType, typename InterType, typename StoreType>
    static void KNNRec(Node* T, const Point& q, const DimsType DIM,
                       kBoundedQueue<Point, StoreType>& bq, const Box& bx,
                       size_t& vis_node_num);

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
