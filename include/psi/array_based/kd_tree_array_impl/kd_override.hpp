#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_

namespace psi {
namespace array_based {

#define KDTREEARRAY_TEMPLATE                                          \
  template <typename Point, typename SplitRule, typename LeafAugType, \
            typename InteriorAugType, uint_fast8_t kSkHeight,         \
            uint_fast8_t kImbaRatio>
#define KDTREEARRAY_CLASS \
  KdTreeArray<Point, SplitRule, LeafAugType, InteriorAugType, kSkHeight, \
              kImbaRatio>

//==============================================================================
// OVERRIDE IMPLEMENTATIONS
//==============================================================================

KDTREEARRAY_TEMPLATE
constexpr void KDTREEARRAY_CLASS::DeleteTree() {
  // Array-based: just clear the vectors
  nodes_.clear();
  leaf_points_.clear();
  this->root_ = BT::NULL_INDEX;
  this->num_points_ = 0;
  this->num_nodes_ = 0;
}

KDTREEARRAY_TEMPLATE
template <typename Range>
void KDTREEARRAY_CLASS::Flatten(Range&& Out) {
  // TODO: Implement
  // Traverse tree and collect all points
  throw std::runtime_error("KdTreeArray::Flatten not yet implemented");
}

KDTREEARRAY_TEMPLATE
auto KDTREEARRAY_CLASS::RangeCount(Box const& query_box) {
  // TODO: Implement
  // Count points in range
  throw std::runtime_error("KdTreeArray::RangeCount not yet implemented");
}

KDTREEARRAY_TEMPLATE
template <typename Range>
auto KDTREEARRAY_CLASS::RangeQuery(Box const& query_box, Range&& Out) {
  // TODO: Implement
  // Query points in range
  throw std::runtime_error("KdTreeArray::RangeQuery not yet implemented");
}

KDTREEARRAY_TEMPLATE
template <typename Range>
auto KDTREEARRAY_CLASS::KNN(NodeIndex idx, Point const& q,
                            kBoundedQueue<Point, Range>& bq) {
  // TODO: Implement
  // K-nearest neighbor query using array indices
  throw std::runtime_error("KdTreeArray::KNN not yet implemented");
}

#undef KDTREEARRAY_TEMPLATE
#undef KDTREEARRAY_CLASS

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_
