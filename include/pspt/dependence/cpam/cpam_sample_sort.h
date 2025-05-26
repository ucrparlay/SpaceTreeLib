// This file is basically the cache-oblivious sorting algorithm from:
//
// Low depth cache-oblivious algorithms.
// Guy E. Blelloch, Phillip B. Gibbons and  Harsha Vardhan Simhadri.
// Proc. ACM symposium on Parallelism in algorithms and architectures (SPAA),
// 2010

#ifndef CPAM_SAMPLE_SORT_H
#define CPAM_SAMPLE_SORT_H

#include <cassert>
#include <concepts>
#include <cstdio>
#include <iterator>
#include <limits>
#include <type_traits>

#include "parlay/delayed_sequence.h"
#include "parlay/internal/bucket_sort.h"
#include "parlay/internal/quicksort.h"
#include "parlay/internal/sequence_ops.h"
#include "parlay/internal/transpose.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/parallel.h"
#include "parlay/relocation.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/utilities.h"

namespace parlay {
namespace internal {
namespace cpam {

// the following parameters can be tuned
constexpr size_t const QUICKSORT_THRESHOLD = 16384;
constexpr size_t const OVER_SAMPLE = 8;

// generates counts in Sc for the number of keys in Sa between consecutive
// values of Sb. Sa and Sb must be sorted
template <typename InIterator, typename PivotIterator, typename CountIterator,
          typename Compare>
void get_bucket_counts(slice<InIterator, InIterator> sA,
                       slice<PivotIterator, PivotIterator> sB,
                       slice<CountIterator, CountIterator> sC, Compare f) {
  using s_size_t = typename std::iterator_traits<CountIterator>::value_type;

  if (sA.size() == 0 || sB.size() == 0) return;
  for (auto& c : sC) c = 0;
  auto itA = sA.begin();
  auto itB = sB.begin();
  auto itC = sC.begin();
  while (true) {
    while (f(*itA, *itB)) {
      assert(itC != sC.end());
      (*itC)++;
      if (++itA == sA.end()) return;
    }
    itB++;
    itC++;
    if (itB == sB.end()) break;
    if (!(f(*(itB - 1), *itB))) {
      while (!f(*itB, *itA)) {
        assert(itC != sC.end());
        (*itC)++;
        if (++itA == sA.end()) return;
      }
      itB++;
      itC++;
      if (itB == sB.end()) break;
    }
  }
  assert(itC != sC.end());
  *itC = static_cast<s_size_t>(sA.end() - itA);
}

template <typename Iterator, typename Compare>
void seq_sort_inplace(slice<Iterator, Iterator> A, Compare const& less,
                      bool stable) {
  using value_type = typename slice<Iterator, Iterator>::value_type;
  if (((sizeof(value_type) > 8) || std::is_pointer<value_type>::value))
    if (!stable)
      quicksort(A.begin(), A.size(), less);
    else
      bucket_sort(A, less, true);  // merge_sort_inplace(A, less);
  else
    bucket_sort(A, less, stable);
}

template <typename filling_curve_t, typename assignment_tag,
          typename InIterator, typename OutIterator, typename Compare>
void seq_sort_(slice<InIterator, InIterator> In,
               slice<OutIterator, OutIterator> Out, Compare const& less,
               bool stable = false) {
  using InputType = typename slice<InIterator, InIterator>::value_type;
  using OutputType = typename slice<OutIterator, OutIterator>::value_type;

  size_t l = In.size();
  for (size_t j = 0; j < l; j++) {
    // assign_dispatch(Out[j], In[j], assignment_tag());
    In[j].SetAugMember(filling_curve_t::Encode(In[j]));

    if constexpr (std::same_as<OutputType, InputType>) {
      // NOTE: output is same as input
      Out[j] = In[j];
    } else if constexpr (std::same_as<OutputType, typename InputType::AT>) {
      // NOTE: output is same as augtype(Key)
      Out[j] = In[j].aug;
    } else if constexpr (std::same_as<typename OutputType::first_type,
                                      typename InputType::AT>) {
      // NOTE: this assumes the output is a pair
      Out[j] = std::make_pair(In[j].aug, &In[j]);
    } else {
      static_assert(false, "Unsupported output type");
    }
  }
  seq_sort_inplace(Out, less, stable);
}

// Copying version of sample sort. This one makes copies of the input
// elements when sorting them into the output. Roughly (\sqrt{n})
// additional copies are also made to copy the pivots. This one can
// be stable.
template <typename filling_curve_t, typename s_size_t = size_t,
          typename InIterator, typename OutIterator, typename Compare>
void sample_sort_(slice<InIterator, InIterator> In,
                  slice<OutIterator, OutIterator> Out, Compare const& less,
                  bool stable = false) {
  // std::cout << "this sample_sort_" << std::endl;
  using value_type = typename slice<InIterator, InIterator>::value_type;
  using output_value_type =
      typename slice<OutIterator, OutIterator>::value_type;

  size_t n = In.size();

  if (n < QUICKSORT_THRESHOLD) {
    seq_sort_<filling_curve_t, uninitialized_copy_tag>(In, Out, less, stable);
  } else {
    // The larger these are, the more comparisons are done but less
    // overhead for the transpose.
    size_t bucket_quotient = 4;
    size_t block_quotient = 4;
    if constexpr (std::is_pointer<value_type>::value) {
      bucket_quotient = 2;
      block_quotient = 3;
    } else if constexpr (sizeof(value_type) > 8) {
      bucket_quotient = 3;
      block_quotient = 3;
    }
    size_t sqrt = static_cast<size_t>(std::sqrt(n));
    size_t num_blocks = size_t{1} << log2_up((sqrt / block_quotient) + 1);
    size_t block_size = ((n - 1) / num_blocks) + 1;
    size_t num_buckets = (sqrt / bucket_quotient) + 1;
    size_t sample_set_size = num_buckets * OVER_SAMPLE;
    size_t m = num_blocks * num_buckets;

    // generate "random" samples with oversampling
    // auto sample_set = sequence<value_type>::from_function(
    //     sample_set_size, [&](size_t i) { return In[hash64(i) % n]; });
    auto sample_set = sequence<
        output_value_type>::from_function(sample_set_size, [&](size_t i) {
      auto pt = &In[hash64(i) % n];
      pt->SetAugMember(filling_curve_t::Encode(*pt));

      if constexpr (std::same_as<output_value_type, value_type>) {
        return *pt;
      } else if constexpr (std::same_as<output_value_type,
                                        typename value_type::AT>) {
        return pt->aug;
      } else if constexpr (std::same_as<typename output_value_type::first_type,
                                        typename value_type::AT>) {
        // NOTE: this assumes the output is a pair
        return std::make_pair(pt->aug, pt);
      } else {
        static_assert(false, "Unsupported output type");
      }
    });
    // sort the samples
    quicksort(sample_set.begin(), sample_set_size, less);

    // subselect samples at even stride
    // auto pivots = sequence<value_type>::from_function(
    //     num_buckets - 1, [&](size_t i) { return sample_set[OVER_SAMPLE * i];
    //     });
    auto pivots = sequence<output_value_type>::from_function(
        num_buckets - 1, [&](size_t i) { return sample_set[OVER_SAMPLE * i]; });

    auto Tmp = uninitialized_sequence<output_value_type>(n);

    // sort each block and merge with samples to get counts for each bucket
    auto counts = sequence<s_size_t>::uninitialized(m + 1);
    counts[m] = 0;
    sliced_for(n, block_size, [&](size_t i, size_t start, size_t end) {
      seq_sort_<filling_curve_t, uninitialized_copy_tag>(
          In.cut(start, end), make_slice(Tmp).cut(start, end), less, stable);
      get_bucket_counts(
          make_slice(Tmp).cut(start, end), make_slice(pivots),
          make_slice(counts).cut(i * num_buckets, (i + 1) * num_buckets), less);
    });

    // move data from blocks to buckets
    auto bucket_offsets = transpose_buckets<uninitialized_relocate_tag>(
        Tmp.begin(), Out.begin(), counts, n, block_size, num_blocks,
        num_buckets);

    // sort within each bucket
    parallel_for(
        0, num_buckets,
        [&](size_t i) {
          size_t start = bucket_offsets[i];
          size_t end = bucket_offsets[i + 1];

          // buckets need not be sorted if two consecutive pivots are equal
          if (i == 0 || i == num_buckets - 1 ||
              less(pivots[i - 1], pivots[i])) {
            seq_sort_inplace(Out.cut(start, end), less, stable);
          }
        },
        1);
  }
}

template <typename filling_curve_t, typename key_entry_pointer,
          typename InputIterator, typename Compare>
auto cpam_sample_sort(slice<InputIterator, InputIterator> A,
                      Compare const& less, bool stable = false) {
  using input_value_type =
      typename slice<InputIterator, InputIterator>::value_type;

  sequence<key_entry_pointer> R =
      sequence<key_entry_pointer>::uninitialized(A.size());
  if (A.size() < (std::numeric_limits<unsigned int>::max)()) {
    sample_sort_<filling_curve_t, unsigned int>(A, make_slice(R), less, stable);
  } else {
    sample_sort_<filling_curve_t, size_t>(A, make_slice(R), less, stable);
  }
  return R;
}

}  // namespace cpam
}  // namespace internal
}  // namespace parlay

#endif  // CPAM_SAMPLE_SORT_H
