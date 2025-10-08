#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_

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
// BUILD OPERATIONS
//==============================================================================

KDTREEARRAY_TEMPLATE
template <typename Range>
void KDTREEARRAY_CLASS::Build(Range&& In) {
  // TODO: Implement
  // 1. Estimate node count and reserve space
  // 2. Build tree recursively
  // 3. Compact arrays if needed
  
  // Placeholder
  throw std::runtime_error("KdTreeArray::Build not yet implemented");
}

KDTREEARRAY_TEMPLATE
void KDTREEARRAY_CLASS::Build_(Slice In) {
  // TODO: Implement
  // Main build implementation
  
  throw std::runtime_error("KdTreeArray::Build_ not yet implemented");
}

KDTREEARRAY_TEMPLATE
auto KDTREEARRAY_CLASS::BuildRecursive(Slice In, Slice Out, DimsType dim,
                                        Box const& bx) -> NodeIndex {
  // TODO: Implement
  // Recursive build using array indices
  // Similar to pointer-based version but returns NodeIndex
  
  return BT::NULL_INDEX;
}

KDTREEARRAY_TEMPLATE
auto KDTREEARRAY_CLASS::SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                                              Box const& bx) -> NodeIndex {
  // TODO: Implement
  // Serial build for small subtrees
  
  return BT::NULL_INDEX;
}

KDTREEARRAY_TEMPLATE
void KDTREEARRAY_CLASS::PickPivots(Slice In, size_t const& n,
                                   SplitterSeq& pivots, DimsType const dim,
                                   BoxSeq& box_seq, Box const& bx) {
  // TODO: Implement
  // Pivot selection logic (can reuse from pointer-based)
}

KDTREEARRAY_TEMPLATE
void KDTREEARRAY_CLASS::DivideRotate(Slice In, SplitterSeq& pivots,
                                     DimsType dim, BucketType idx,
                                     BoxSeq& box_seq, Box const& bx) {
  // TODO: Implement
  // Partitioning logic (can reuse from pointer-based)
}

//==============================================================================
// NODE ALLOCATION
//==============================================================================

KDTREEARRAY_TEMPLATE
auto KDTREEARRAY_CLASS::AllocNode() -> NodeIndex {
  nodes_.emplace_back();
  return nodes_.size() - 1;
}

KDTREEARRAY_TEMPLATE
auto KDTREEARRAY_CLASS::AllocLeafNode(Slice In) -> NodeIndex {
  NodeIndex idx = AllocNode();
  CompactNode& node = nodes_[idx];
  
  node.is_leaf = true;
  node.size = In.size();
  node.leaf_start_idx = leaf_points_.size();
  node.leaf_count = In.size();
  
  // Copy points to contiguous storage
  for (auto const& p : In) {
    leaf_points_.push_back(p);
  }
  
  return idx;
}

KDTREEARRAY_TEMPLATE
void KDTREEARRAY_CLASS::FreeNode(NodeIndex idx) {
  // For array-based, we don't actually free individual nodes
  // Freeing happens during tree rebuild/compaction
}

#undef KDTREEARRAY_TEMPLATE
#undef KDTREEARRAY_CLASS

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BUILD_TREE_HPP_
