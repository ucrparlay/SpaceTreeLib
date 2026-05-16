#include <cstdint>
#include <iostream>
#include <vector>

#include "partition/vqsort_partition.hpp"

int main() {
  std::vector<uint64_t> values = {
      19, 4, 7, 31, 12, 1, 9, 22, 16, 5, 27, 8, 14, 3, 25, 11,
  };
  const uint64_t pivot = 12;

  const size_t bound =
      psi::partition::PartitionAscending(values.data(), values.size(), pivot);

  std::cout << "partition bound = " << bound << '\n';
  std::cout << "[<= pivot]:";
  for (size_t i = 0; i < bound; ++i) std::cout << ' ' << values[i];
  std::cout << "\n[> pivot]:";
  for (size_t i = bound; i < values.size(); ++i) std::cout << ' ' << values[i];
  std::cout << '\n';

  std::cout << "uint64_t scalar fallback cutoff = "
            << psi::partition::VQSortPartitionBaseCaseKeysAscending<uint64_t>()
            << " keys\n";
}
