#pragma once

#include <algorithm>
#include <cstddef>
#include <type_traits>

// This wrapper intentionally uses vqsort's internal partition kernel. Defining
// VQSORT_ONLY_STATIC avoids requiring the dynamic-dispatch vqsort object files.
#ifndef VQSORT_ONLY_STATIC
#define VQSORT_ONLY_STATIC 1
#endif

#include "hwy/base.h"
#include "hwy/contrib/sort/order.h"
#include "hwy/contrib/sort/vqsort-inl.h"

namespace psi::partition {

namespace detail {

template <typename T>
inline bool IsRightPartition(const T& value, const T& pivot,
                             hwy::SortAscending) {
  return pivot < value;
}

template <typename T>
inline bool IsRightPartition(const T& value, const T& pivot,
                             hwy::SortDescending) {
  return value < pivot;
}

// Same boundary semantics as vqsort's two-way partition:
// ascending  => [<= pivot][> pivot]
// descending => [>= pivot][< pivot]
template <typename T, class Order>
size_t ScalarPartition(T* keys, size_t num, const T& pivot, Order order) {
  size_t left = 0;
  size_t right = num;

  while (left < right) {
    while (left < right && !IsRightPartition(keys[left], pivot, order)) {
      ++left;
    }
    while (left < right && IsRightPartition(keys[right - 1], pivot, order)) {
      --right;
    }
    if (left < right) {
      --right;
      std::swap(keys[left], keys[right]);
      ++left;
    }
  }

  return left;
}

template <class Traits>
constexpr size_t LanesPerKey() {
  return Traits().LanesPerKey();
}

}  // namespace detail

template <typename Key, class Order>
size_t VQSortPartition(Key* HWY_RESTRICT keys, size_t num_keys,
                       const Key& pivot, Order order) {
  static_assert(std::is_same_v<Order, hwy::SortAscending> ||
                    std::is_same_v<Order, hwy::SortDescending>,
                "Order must be hwy::SortAscending or hwy::SortDescending");

  namespace hn = hwy::HWY_NAMESPACE;

  const hn::detail::MakeTraits<Key, Order> st;
  using Traits = decltype(st);
  using LaneType = typename Traits::LaneType;
  constexpr size_t kLanesPerKey = detail::LanesPerKey<Traits>();

  const hn::SortTag<LaneType> d;
  const size_t lanes = hn::Lanes(d);
  const size_t num_lanes = num_keys * kLanesPerKey;
  const size_t base_case_num =
      hwy::SortConstants::BaseCaseNumLanes<kLanesPerKey>(lanes);

  if (HWY_UNLIKELY(num_lanes <= base_case_num)) {
    return detail::ScalarPartition(keys, num_keys, pivot, order);
  }

  HWY_ALIGN LaneType pivot_lanes[kLanesPerKey];
  hwy::CopyBytes(&pivot, pivot_lanes, sizeof(Key));
  const auto pivot_vec = st.SetKey(d, pivot_lanes);

  HWY_ALIGN LaneType buf[hwy::SortConstants::BufBytes<
                             LaneType, kLanesPerKey>(HWY_MAX_BYTES) /
                         sizeof(LaneType)];

  const size_t bound_lanes = hn::detail::Partition(
      d, st, reinterpret_cast<LaneType* HWY_RESTRICT>(keys), num_lanes,
      pivot_vec, buf);

  return bound_lanes / kLanesPerKey;
}

template <typename Key>
size_t PartitionAscending(Key* HWY_RESTRICT keys, size_t num_keys,
                          const Key& pivot) {
  return VQSortPartition(keys, num_keys, pivot, hwy::SortAscending());
}

template <typename Key>
size_t PartitionDescending(Key* HWY_RESTRICT keys, size_t num_keys,
                           const Key& pivot) {
  return VQSortPartition(keys, num_keys, pivot, hwy::SortDescending());
}

// Exposed so callers can choose their own small-range path before calling the
// SIMD wrapper. The value is in keys, not Highway lanes.
template <typename Key, class Order>
size_t VQSortPartitionBaseCaseKeys(Order) {
  namespace hn = hwy::HWY_NAMESPACE;

  const hn::detail::MakeTraits<Key, Order> st;
  using Traits = decltype(st);
  using LaneType = typename Traits::LaneType;
  constexpr size_t kLanesPerKey = detail::LanesPerKey<Traits>();

  const hn::SortTag<LaneType> d;
  const size_t lanes = hn::Lanes(d);
  return hwy::SortConstants::BaseCaseNumLanes<kLanesPerKey>(lanes) /
         kLanesPerKey;
}

template <typename Key>
size_t VQSortPartitionBaseCaseKeysAscending() {
  return VQSortPartitionBaseCaseKeys<Key>(hwy::SortAscending());
}

template <typename Key>
size_t VQSortPartitionBaseCaseKeysDescending() {
  return VQSortPartitionBaseCaseKeys<Key>(hwy::SortDescending());
}

}  // namespace psi::partition
