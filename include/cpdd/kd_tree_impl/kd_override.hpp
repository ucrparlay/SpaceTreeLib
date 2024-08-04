#include <type_traits>
#include <utility>
#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::KNN(Node* T, const Point& q,
                                         const DimsType DIM,
                                         kBoundedQueue<Point, Range>& bq,
                                         const Box& bx, size_t& vis_node_num) {
    BT::template KNNBinary<Leaf, Interior>(T, q, DIM, bq, bx, vis_node_num);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::Flatten(Range&& Out) {
    BT::template FlattenRec<Leaf, Interior>(this->root_,
                                            parlay::make_slice(Out));
}

template<typename Point, typename SplitRule, uint8_t kBDO>
size_t KdTree<Point, SplitRule, kBDO>::RangeCount(const Box& bx) {
    size_t vis_leaf_num = 0, vis_inter_num = 0;
    return BT::template RangeCountRectangle<Leaf, Interior>(
        this->root_, bx, this->tree_box_, vis_leaf_num, vis_inter_num);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
size_t KdTree<Point, SplitRule, kBDO>::RangeCount(const Circle& cl) {
    return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
                                                         this->tree_box_);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
template<typename Range>
size_t KdTree<Point, SplitRule, kBDO>::RangeQuery(const Box& query_box,
                                                  Range&& Out) {
    size_t s = 0;
    BT::template RangeQuerySerialRecursive<Leaf, Interior>(
        this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_);
    return s;
}

template<typename Point, typename SplitRule, uint8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::DeleteTree() {
    this->template DeleteTreeWrapper<Leaf, Interior, BinaryInteriorTag>();
}

}  // namespace cpdd
