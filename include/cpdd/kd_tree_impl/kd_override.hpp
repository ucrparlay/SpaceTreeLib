#include <type_traits>
#include <utility>
#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kBDO>
template<typename Range>
auto KdTree<Point, SplitRule, kBDO>::KNN(Node* T, const Point& q,
                                         kBoundedQueue<Point, Range>& bq) {
    KNNLogger logger;
    BT::template KNNBinary<Leaf, Interior>(T, q, bq, this->tree_box_, logger);
    return logger;
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::Flatten(Range&& Out) {
    BT::template FlattenRec<Leaf, Interior>(this->root_,
                                            parlay::make_slice(Out));
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
size_t KdTree<Point, SplitRule, kBDO>::RangeCount(const Box& bx) {
    return BT::template RangeCountRectangle<Leaf, Interior>(this->root_, bx,
                                                            this->tree_box_);
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
size_t KdTree<Point, SplitRule, kBDO>::RangeCount(const Circle& cl) {
    return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
                                                         this->tree_box_);
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
template<typename Range>
size_t KdTree<Point, SplitRule, kBDO>::RangeQuery(const Box& query_box,
                                                  Range&& Out) {
    size_t s = 0;
    BT::template RangeQuerySerialRecursive<Leaf, Interior>(
        this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_);
    return s;
}

template<typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::DeleteTree() {
    BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace cpdd
