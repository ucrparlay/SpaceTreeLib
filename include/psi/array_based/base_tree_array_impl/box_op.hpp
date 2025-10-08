#ifndef PSI_ARRAY_BASED_BASE_TREE_ARRAY_IMPL_BOX_OP_HPP_
#define PSI_ARRAY_BASED_BASE_TREE_ARRAY_IMPL_BOX_OP_HPP_

namespace psi {
namespace array_based {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto BaseTreeArray<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(
    Range const& In) -> Box {
  if (In.size() == 0) {
    return GetEmptyBox();
  }

  auto min_coords = BasicPoint::GetMax();
  auto max_coords = BasicPoint::GetMin();

  for (auto const& p : In) {
    for (DimsType d = 0; d < kDim; ++d) {
      min_coords.pnt[d] = std::min(min_coords.pnt[d], p.pnt[d]);
      max_coords.pnt[d] = std::max(max_coords.pnt[d], p.pnt[d]);
    }
  }

  return Box(min_coords, max_coords);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
auto BaseTreeArray<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(
    Box const& a, Box const& b) -> Box {
  auto min_coords = a.first;
  auto max_coords = a.second;

  for (DimsType d = 0; d < kDim; ++d) {
    min_coords.pnt[d] = std::min(min_coords.pnt[d], b.first.pnt[d]);
    max_coords.pnt[d] = std::max(max_coords.pnt[d], b.second.pnt[d]);
  }

  return Box(min_coords, max_coords);
}

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_BASE_TREE_ARRAY_IMPL_BOX_OP_HPP_
