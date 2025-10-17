#ifndef PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_
#define PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_

#include "array_view/kd_tree_array.h"
namespace psi {
namespace array_view {

//==============================================================================
// OVERRIDE IMPLEMENTATIONS
//==============================================================================

template <typename TypeTrait>
constexpr void KdTreeArray<TypeTrait>::DeleteTree() {
  // BUG: needs to handle the DeleteTree
  // this->inner_tree_seq_.clear();
  // this->leaf_seq_.clear();
  return;
}

template <typename TypeTrait>
template <typename Range>
void KdTreeArray<TypeTrait>::Flatten(Range&& Out) {
  // TODO: Implement
  // Traverse tree and collect all points
  throw std::runtime_error("KdTreeArray::Flatten not yet implemented");
}

template <typename TypeTrait>
auto KdTreeArray<TypeTrait>::RangeCount(Box const& query_box) {
  // TODO: Implement
  // Count points in range
  RangeQueryLogger logger;
  throw std::runtime_error("KdTreeArray::RangeCount not yet implemented");
  size_t size = 0;
  return std::make_pair(size, logger);
}

template <typename TypeTrait>
template <typename Range>
auto KdTreeArray<TypeTrait>::RangeQuery(Box const& query_box, Range&& Out) {
  // TODO: Implement
  // Query points in range
  RangeQueryLogger logger;
  throw std::runtime_error("KdTreeArray::RangeQuery not yet implemented");
  size_t size = 0;
  return std::make_pair(size, logger);
}

template <typename TypeTrait>
template <typename Range>
auto KdTreeArray<TypeTrait>::KNN(NodeIndex idx, Point const& q,
                                 kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  throw std::runtime_error("KdTreeArray::KNN not yet implemented");
  return logger;
}

}  // namespace array_view
}  // namespace psi

#endif  // PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_OVERRIDE_HPP_
