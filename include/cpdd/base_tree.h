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
    using slice = parlay::slice<Point*, Point*>;
    using Points = parlay::sequence<Point>;
    using PointsIter = typename parlay::sequence<Point>::iterator;
    using Splitter = std::pair<Coord, DimsType>;
    using SplitterSeq = parlay::sequence<Splitter>;
    using Box = std::pair<Point, Point>;
    using BoxSeq = parlay::sequence<Box>;
    using Circle = std::pair<Point, Coord>;

    using NodeBox = std::pair<node*, Box>;
    using NodeTag = std::pair<node*, uint_fast8_t>;
    using NodeTagSeq = parlay::sequence<NodeTag>;
    using TagNodes = parlay::sequence<BallsType>;  //*index by tag
    //@ Const variables
    //@ uint32t handle up to 4e9 at least
    //! bucket num should smaller than 1<<8 to handle type overflow

    // TODO wrap the variables using std::getenv()
    static constexpr BucketType kBuildDepthOnce = 6;  //* last layer is leaf
    static constexpr BucketType kPivotNum = (1 << kBuildDepthOnce) - 1;
    static constexpr BucketType kBucketNum = 1 << kBuildDepthOnce;
    //@ tree structure
    static constexpr uint_fast8_t kLeaveWrap = 32;
    static constexpr uint_fast8_t kThinLeaveWrap = 24;
    static constexpr uint_fast16_t kSerialBuildCutoff = 1 << 10;
    //@ block param in Partition
    static constexpr uint_fast8_t kLog2Base = 10;
    static constexpr uint_fast16_t kBlockSize = 1 << kLog2Base;
    //@ reconstruct weight threshold
    static constexpr uint_fast8_t kInbalanceRatio = 30;

    inline size_t GetImbalanceRatio();
    inline bool ImbalanceNode(const size_t l, const size_t n);

    enum SplitRule { kMaxStretchDim, kRotateDim };

    //@ array based inner tree for batch insertion and deletion
    struct InnerTree;

    //@ Box operations
    static inline bool LegalBox(const Box& bx);
    static inline bool WithinBox(const Box& a, const Box& b);
    static inline bool WithinBox(const Point& p, const Box& bx);
    static inline bool BoxIntersectBox(const Box& a, const Box& b);
    static inline Box GetEmptyBox();
    static Box GetBox(const Box& x, const Box& y);
    static Box GetBox(slice V);
    static Box GetBox(node* T);

    static inline bool WithinCircle(const Box& bx, const Circle& cl);
    static inline bool WithinCircle(const Point& p, const Circle& cl);
    static inline bool CircleIntersectBox(const Circle& cl, const Box& bx);

    //@ dimensionality
    inline DimsType PickRebuildDim(const node* T, const DimsType d, const DimsType DIM);
    static inline DimsType PickMaxStretchDim(const Box& bx, const DimsType DIM);

    //@ Parallel KD tree cores
    //@ Build
    void DivideRotate(slice In, SplitterSeq& pivots, DimsType dim, BucketType idx, BucketType deep, BucketType& bucket,
                      const DimsType DIM, BoxSeq& boxs, const Box& bx);
    void PickPivots(slice In, const size_t& n, SplitterSeq& pivots, const DimsType dim, const DimsType DIM,
                    BoxSeq& boxs, const Box& bx);
    static inline BucketType FindBucket(const Point& p, const SplitterSeq& pivots);
    static void Partition(slice A, slice B, const size_t n, const SplitterSeq& pivots,
                          parlay::sequence<BallsType>& sums);
    static node* BuildInnerTree(BucketType idx, SplitterSeq& pivots, parlay::sequence<node*>& treeNodes);
    PointsIter SerialPartition(slice In, DimsType d);

    virtual void Build(slice In, const DimsType DIM) = 0;
    // virtual node* serial_build_recursive(slice In, slice Out, DimsType dim,
    //                                      const DimsType DIM, const Box& bx) =
    //                                      0;
    // virtual node* build_recursive(slice In, slice Out, DimsType dim,
    //                               const DimsType DIM, const Box& bx) = 0;

    //@ random support
    static inline uint64_t Hash64(uint64_t u);

    //@ batch helpers:
    template<typename Slice>
    static void Flatten(node* T, Slice Out, bool granularity = true);

    void FlattenAndDelete(node* T, slice Out);
    static void SeievePoints(slice A, slice B, const size_t n, const NodeTagSeq& tags,
                             parlay::sequence<BallsType>& sums, const BucketType tagsNum);
    static inline BucketType RetriveTag(const Point& p, const NodeTagSeq& tags);
    static node* UpdateInnerTree(BucketType idx, const NodeTagSeq& tags, parlay::sequence<node*>& treeNodes,
                                 BucketType& p, const TagNodes& rev_tag);

    virtual void DeleteTree() = 0;

    template<typename leaf_node_type, typename interior_node_type>
    void DeleteTreeWrapper();

    template<typename leaf_node_type, typename interior_node_type>
    static void DeleteTreeRecursive(node* T, bool granularity = true);

    //@ batch insert
    node* RebuildWithInsert(node* T, slice In, const DimsType d, const DimsType DIM);
    static inline void UpdateInterior(node* T, node* L, node* R);
    void BatchInsert(slice In, const DimsType DIM);
    node* BatchInsertRecursive(node* T, slice In, slice Out, DimsType d, const DimsType DIM);

    //@ batch delete
    NodeBox RebuildAfterDelete(node* T, const DimsType d, const DimsType DIM);
    void BatchDelete(slice In, const DimsType DIM);
    NodeBox BatchDeleteRecursive(node* T, slice In, slice Out, DimsType d, const DimsType DIM, bool hasTomb);
    NodeBox DeleteInnerTree(BucketType idx, const NodeTagSeq& tags, parlay::sequence<NodeBox>& treeNodes, BucketType& p,
                            const TagNodes& rev_tag, const DimsType d, const DimsType DIM);

    //@ validations
    template<typename interior>
    bool CheckBox(node* T, const Box& bx);

    template<typename interior>
    size_t CheckSize(node* T);

    template<typename interior>
    void CheckTreeSameSequential(node* T, int dim, const int& DIM);

    template<typename interior>
    void Validate(const DimsType DIM);

    template<typename interior>
    size_t GetTreeHeight();

    template<typename interior>
    size_t GetMaxTreeDepth(node* T, size_t deep);

    template<typename interior>
    double GetAveTreeHeight();

    template<typename interior>
    size_t CountTreeNodesNum(node* T);

    template<typename interior>
    void CountTreeHeights(node* T, size_t deep, size_t& idx, parlay::sequence<size_t>& heights);

    //@ kdtree interfaces
    inline void SetRoot(node* _root) { this->root_ = _root; }

    inline node* GetRoot() { return this->root_; }

    inline Box GetRootBox() { return this->tree_box_; }

   protected:
    node* root_ = nullptr;
    parlay::internal::timer timer;
    // SplitRule split_rule_ = kRotateDim;
    SplitRule split_rule_ = kMaxStretchDim;
    Box tree_box_;
};

}  // namespace cpdd

#include "base_tree_impl/box_op.hpp"
// #include "base_tree_impl/validation.hpp"
#include "base_tree_impl/delete_tree.hpp"
#include "base_tree_impl/dimensinality.hpp"
#include "base_tree_impl/random_support.hpp"
