#include <type_traits>
#include <utility>
#include "../kd_tree.h"

namespace cpdd {
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
