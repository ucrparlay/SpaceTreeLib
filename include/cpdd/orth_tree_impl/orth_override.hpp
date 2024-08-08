#include <tuple>
#include <type_traits>
#include <utility>
#include "../orth_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
template<typename Range>
auto OrthTree<Point, SplitRule, kMD, kBDO>::KNN(
    Node* T, const Point& q, const DimsType DIM,
    kBoundedQueue<Point, Range>& bq) {
    size_t vis_node_num = 0, generate_box_num = 0, check_box_num = 0;
    BT::template KNNBinary<Leaf, KdInteriorNode>(T, q, DIM, bq, this->tree_box_,
                                                 vis_node_num, generate_box_num,
                                                 check_box_num);
    // BT::template KNNMulti<Leaf, Interior>(T, q, DIM, bq, bx, vis_node_num);
    return std::make_tuple(vis_node_num, generate_box_num, check_box_num);
}

template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
template<typename Range>
void OrthTree<Point, SplitRule, kMD, kBDO>::Flatten(Range&& Out) {
    BT::template FlattenRec<Leaf, Interior>(this->root_,
                                            parlay::make_slice(Out));
}

template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
size_t OrthTree<Point, SplitRule, kMD, kBDO>::RangeCount(const Box& bx) {
    size_t vis_leaf_num = 0, vis_inter_num = 0;
    return BT::template RangeCountRectangle<Leaf, Interior>(
        this->root_, bx, this->tree_box_, vis_leaf_num, vis_inter_num);
}

template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
size_t OrthTree<Point, SplitRule, kMD, kBDO>::RangeCount(const Circle& cl) {
    return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
                                                         this->tree_box_);
}

template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
template<typename Range>
size_t OrthTree<Point, SplitRule, kMD, kBDO>::RangeQuery(const Box& query_box,
                                                         Range&& Out) {
    size_t s = 0;
    BT::template RangeQuerySerialRecursive<Leaf, Interior>(
        this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_);
    return s;
}

template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::DeleteTree() {
    BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace cpdd
