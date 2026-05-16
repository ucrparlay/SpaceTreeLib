#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "hwy/base.h"
#include "partition/vqsort_partition.hpp"

namespace {

template <typename T>
bool FullLess(const T& a, const T& b) {
  return a < b;
}

template <>
bool FullLess(const hwy::K32V32& a, const hwy::K32V32& b) {
  return a.key == b.key ? a.value < b.value : a.key < b.key;
}

template <>
bool FullLess(const hwy::K64V64& a, const hwy::K64V64& b) {
  return a.key == b.key ? a.value < b.value : a.key < b.key;
}

template <>
bool FullLess(const hwy::uint128_t& a, const hwy::uint128_t& b) {
  return a.hi == b.hi ? a.lo < b.lo : a.hi < b.hi;
}

template <typename T>
bool IsRightAscending(const T& value, const T& pivot) {
  return pivot < value;
}

template <typename T>
bool IsRightDescending(const T& value, const T& pivot) {
  return value < pivot;
}

template <typename T>
bool SameMultiset(std::vector<T> before, std::vector<T> after) {
  auto less = [](const T& a, const T& b) { return FullLess(a, b); };
  std::sort(before.begin(), before.end(), less);
  std::sort(after.begin(), after.end(), less);
  return std::equal(before.begin(), before.end(), after.begin(),
                    [](const T& a, const T& b) {
                      return !FullLess(a, b) && !FullLess(b, a);
                    });
}

template <typename T, typename IsRight>
bool VerifyPartition(const std::vector<T>& before, const std::vector<T>& after,
                     size_t bound, const T& pivot, IsRight is_right) {
  if (before.size() != after.size() || bound > after.size()) return false;
  if (!SameMultiset(before, after)) return false;

  for (size_t i = 0; i < bound; ++i) {
    if (is_right(after[i], pivot)) return false;
  }
  for (size_t i = bound; i < after.size(); ++i) {
    if (!is_right(after[i], pivot)) return false;
  }
  return true;
}

uint64_t Mix(uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}

template <typename T>
T MakeValue(size_t i);

template <>
uint32_t MakeValue<uint32_t>(size_t i) {
  return static_cast<uint32_t>(Mix(i) % 1024);
}

template <>
int64_t MakeValue<int64_t>(size_t i) {
  return static_cast<int64_t>(Mix(i) % 2048) - 1024;
}

template <>
hwy::uint128_t MakeValue<hwy::uint128_t>(size_t i) {
  return hwy::uint128_t{Mix(i * 3 + 1), Mix(i * 3 + 2) % 1024};
}

template <>
hwy::K32V32 MakeValue<hwy::K32V32>(size_t i) {
  return hwy::K32V32{static_cast<uint32_t>(i),
                     static_cast<uint32_t>(Mix(i) % 1024)};
}

template <>
hwy::K64V64 MakeValue<hwy::K64V64>(size_t i) {
  return hwy::K64V64{Mix(i * 5 + 1), Mix(i * 5 + 2) % 1024};
}

template <typename T>
std::vector<T> MakeInput(size_t n) {
  std::vector<T> values(n);
  for (size_t i = 0; i < n; ++i) {
    values[i] = MakeValue<T>(i);
  }
  return values;
}

template <typename T>
bool RunOneType(const std::string& name, const T& pivot) {
  bool ok = true;

  for (size_t n : {0ul, 1ul, 7ul, 31ul, 257ul, 4097ul}) {
    auto input = MakeInput<T>(n);
    auto asc = input;
    const size_t asc_bound =
        psi::partition::PartitionAscending(asc.data(), asc.size(), pivot);
    ok &= VerifyPartition(input, asc, asc_bound, pivot, IsRightAscending<T>);

    auto desc = input;
    const size_t desc_bound =
        psi::partition::PartitionDescending(desc.data(), desc.size(), pivot);
    ok &= VerifyPartition(input, desc, desc_bound, pivot, IsRightDescending<T>);
  }

  std::cout << name << ": " << (ok ? "ok" : "failed") << '\n';
  return ok;
}

bool RunAdjacentRanges() {
  constexpr uint64_t kGuard = 0xdeadbeefcafebabeull;
  const size_t n1 = 4097;
  const size_t n2 = 5003;
  const size_t start1 = 17;
  const size_t start2 = start1 + n1;
  std::vector<uint64_t> data(start2 + n2 + 19, kGuard);

  for (size_t i = 0; i < n1; ++i) data[start1 + i] = MakeValue<uint32_t>(i);
  for (size_t i = 0; i < n2; ++i) data[start2 + i] = MakeValue<uint32_t>(i + 99);

  const auto before1 =
      std::vector<uint64_t>(data.begin() + start1, data.begin() + start2);
  const auto before2 =
      std::vector<uint64_t>(data.begin() + start2, data.begin() + start2 + n2);

  const uint64_t pivot1 = 500;
  const size_t bound1 =
      psi::partition::PartitionAscending(data.data() + start1, n1, pivot1);
  const auto after1 =
      std::vector<uint64_t>(data.begin() + start1, data.begin() + start2);
  const auto untouched2 =
      std::vector<uint64_t>(data.begin() + start2, data.begin() + start2 + n2);

  bool ok = VerifyPartition(before1, after1, bound1, pivot1,
                            IsRightAscending<uint64_t>);
  ok &= SameMultiset(before2, untouched2);

  const uint64_t pivot2 = 700;
  const size_t bound2 =
      psi::partition::PartitionAscending(data.data() + start2, n2, pivot2);
  const auto after2 =
      std::vector<uint64_t>(data.begin() + start2, data.begin() + start2 + n2);
  ok &= VerifyPartition(before2, after2, bound2, pivot2,
                        IsRightAscending<uint64_t>);

  for (size_t i = 0; i < start1; ++i) ok &= data[i] == kGuard;
  for (size_t i = start2 + n2; i < data.size(); ++i) ok &= data[i] == kGuard;

  std::cout << "adjacent ranges: " << (ok ? "ok" : "failed") << '\n';
  return ok;
}

}  // namespace

int main() {
  bool ok = true;
  ok &= RunOneType<uint32_t>("uint32_t", uint32_t{512});
  ok &= RunOneType<int64_t>("int64_t", int64_t{0});
  ok &= RunOneType<hwy::K32V32>("K32V32", hwy::K32V32{0, 512});
  ok &= RunOneType<hwy::uint128_t>("uint128_t", hwy::uint128_t{0, 512});
  ok &= RunOneType<hwy::K64V64>("K64V64", hwy::K64V64{0, 512});
  ok &= RunAdjacentRanges();

  const size_t small =
      psi::partition::VQSortPartitionBaseCaseKeysAscending<uint64_t>();
  std::cout << "uint64_t SIMD cutoff: n <= " << small
            << " uses scalar fallback\n";

  return ok ? 0 : 1;
}
