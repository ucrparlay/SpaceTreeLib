#include <type_traits>
#include <utility>
#include "../kd_tree.h"

namespace cpdd {
template<typename Point, typename SplitRule>
template<typename StoreType>
void KdTree<Point, SplitRule>::KNN(Node* T, const Point& q, const DimsType DIM,
                                   kBoundedQueue<Point, StoreType>& bq,
                                   const Box& bx, size_t& vis_node_num) {
    BT::template KNNRec<Leaf, Interior>(T, q, DIM, bq, bx, vis_node_num);
}

template<typename Point, typename SplitRule>
template<typename Range>
void KdTree<Point, SplitRule>::Flatten(Range&& Out) {
    BT::template FlattenRec<Leaf, Interior>(
        this->root_, parlay::make_slice(Out),
        [&](const Interior* TI) { return TI->aug; }, false);
}

template<typename Point, typename SplitRule>
void KdTree<Point, SplitRule>::DeleteTree() {
    this->template DeleteTreeWrapper<
        LeafNode<Point, Slice, BT::kLeaveWrap, std::true_type>,
        InteriorNode<Point, Splitter, bool>>();
}

}  // namespace cpdd
