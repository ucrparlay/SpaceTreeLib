#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_

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
// BATCH DELETE OPERATIONS
//==============================================================================

KDTREEARRAY_TEMPLATE
template <typename Range>
void KDTREEARRAY_CLASS::BatchDelete(Range&& In) {
  // TODO: Implement
  // Strategy: Rebuild tree without deleted points
  throw std::runtime_error("KdTreeArray::BatchDelete not yet implemented");
}

#undef KDTREEARRAY_TEMPLATE
#undef KDTREEARRAY_CLASS

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_
