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
    using Splitter = std::pair<Coord, DimsType>;
    using SplitterSeq = parlay::sequence<Splitter>;
    using Box = std::pair<Point, Point>;
    using BoxSeq = parlay::sequence<Box>;
    using Circle = std::pair<Point, Coord>;

    using NodeBox = std::pair<Node*, Box>;
    using NodeTag = std::pair<Node*, uint_fast8_t>;
    using NodeTagSeq = parlay::sequence<NodeTag>;
    using TagNodes = parlay::sequence<BallsType>;  //*index by tag
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

    // NOTE: split rule
    struct MaxStretchDimTag {};
    struct RotateDimTag {};

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
    inline DimsType PickRebuildDim_(const Node* T, const DimsType d,
                                    const DimsType DIM, MaxStretchDimTag);
    inline DimsType PickRebuildDim_(const Node* T, const DimsType d,
                                    const DimsType DIM, RotateDimTag);
    static inline DimsType PickMaxStretchDim(const Box& bx, const DimsType DIM);

    // NOTE: build tree
    void DivideRotate(Slice In, SplitterSeq& pivots, DimsType dim,
                      BucketType idx, BucketType deep, BucketType& bucket,
                      const DimsType DIM, BoxSeq& boxs, const Box& bx);

    void PickPivots(Slice In, const size_t& n, SplitterSeq& pivots,
                    const DimsType dim, const DimsType DIM, BoxSeq& boxs,
                    const Box& bx);

    static inline BucketType FindBucket(const Point& p,
                                        const SplitterSeq& pivots);

    static void Partition(Slice A, Slice B, const size_t n,
                          const SplitterSeq& pivots,
                          parlay::sequence<BallsType>& sums);

    static Node* BuildInnerTree(BucketType idx, SplitterSeq& pivots,
                                parlay::sequence<Node*>& tree_nodes);

    PointsIter SerialPartition(Slice In, DimsType d);

    virtual void Build_(Slice In, const DimsType DIM) = 0;

    inline uint64_t Hash64(uint64_t u);
    virtual void DeleteTree() = 0;

    template<typename leaf_node_type, typename interior_node_type>
    void DeleteTreeWrapper();

    template<typename leaf_node_type, typename interior_node_type>
    static void DeleteTreeRecursive(Node* T, bool granularity = true);

    //@ validations
    template<typename interior>
    bool CheckBox(Node* T, const Box& bx);

    template<typename interior>
    size_t CheckSize(Node* T);

    template<typename interior>
    void CheckTreeSameSequential(Node* T, int dim, const int& DIM);

    template<typename interior>
    void Validate(const DimsType DIM);

    template<typename interior>
    size_t GetTreeHeight();

    template<typename interior>
    size_t GetMaxTreeDepth(Node* T, size_t deep);

    template<typename interior>
    double GetAveTreeHeight();

    template<typename interior>
    size_t CountTreeNodesNum(Node* T);

    template<typename interior>
    void CountTreeHeights(Node* T, size_t deep, size_t& idx,
                          parlay::sequence<size_t>& heights);

    //@ kdtree interfaces
    inline void SetRoot(Node* _root) { this->root_ = _root; }

    inline Node* GetRoot() { return this->root_; }

    inline Box GetRootBox() { return this->tree_box_; }

   protected:
    Node* root_ = nullptr;
    parlay::internal::timer timer;
    // SplitRule split_rule_ = kRotateDim;
    Box tree_box_;
};

}  // namespace cpdd

#include "base_tree_impl/box_op.hpp"
// #include "base_tree_impl/validation.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/random_support.hpp"
