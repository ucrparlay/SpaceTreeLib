#ifndef PSI_BASE_TREE_IMPL_INNER_TREE_COMMON_HPP_
#define PSI_BASE_TREE_IMPL_INNER_TREE_COMMON_HPP_

#include <cstdint>

#include "base_tree.h"

namespace psi {
namespace inner_tree_detail {

static constexpr int_fast32_t kInvalidBucket = -1;

enum class NodeTagValue : uint_fast8_t {
  kNormal = 0,
  kLeaf = 1,
  kAncestorRebuilt = 2,
  kNeedsRebuild = 3
};

template <typename BucketType>
inline BucketType GetTagValue(NodeTagValue tag, BucketType kBucketNum) {
  return kBucketNum + static_cast<BucketType>(tag);
}

}  // namespace inner_tree_detail
}  // namespace psi

#endif  // PSI_BASE_TREE_IMPL_INNER_TREE_COMMON_HPP_
