// This file is basically the cache-oblivious sorting algorithm from:
//
// Low depth cache-oblivious algorithms.
// Guy E. Blelloch, Phillip B. Gibbons and  Harsha Vardhan Simhadri.
// Proc. ACM symposium on Parallelism in algorithms and architectures (SPAA),
// 2010

#ifndef CPAM_SAMPLE_SORT_H
#define CPAM_SAMPLE_SORT_H

#include <cassert>
#include <chrono>
#include <concepts>
#include <cstdio>
#include <cstdint>
#include <iterator>
#include <limits>
#include <unistd.h>
#include <utility>
#include <type_traits>

#include "parlay/delayed_sequence.h"
#include "parlay/internal/bucket_sort.h"
#include "parlay/internal/quicksort.h"
#include "parlay/internal/sequence_ops.h"
#include "parlay/internal/transpose.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/relocation.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/utilities.h"

#if __has_include("SIMD-Sample-Sort/src/modules/utils/KVPair.hpp")
#include "SIMD-Sample-Sort/src/modules/utils/KVPair.hpp"
#define PSI_CPAM_HAS_KVPAIR 1
#else
#define PSI_CPAM_HAS_KVPAIR 0
#endif

#ifndef PSI_USE_SIMD_SAMPLE_SORT
#define PSI_USE_SIMD_SAMPLE_SORT 0
#endif

#ifndef PSI_CPAM_PRECOMPUTE_CODE_SORT
#define PSI_CPAM_PRECOMPUTE_CODE_SORT 0
#endif

#ifndef PSI_PRINT_CPAM_BUILD_TIMING
#define PSI_PRINT_CPAM_BUILD_TIMING 0
#endif

#if PSI_USE_SIMD_SAMPLE_SORT
#if PSI_CPAM_HAS_KVPAIR && \
    __has_include("SIMD-Sample-Sort/src/two_pass/two_pass_simd.hpp")
#include "SIMD-Sample-Sort/src/two_pass/two_pass_simd.hpp"
#define PSI_CPAM_HAS_SIMD_SAMPLE_SORT 1
#else
#define PSI_CPAM_HAS_SIMD_SAMPLE_SORT 0
#endif
#else
#define PSI_CPAM_HAS_SIMD_SAMPLE_SORT 0
#endif

namespace parlay {
namespace internal {
namespace cpam {

#if PSI_PRINT_CPAM_BUILD_TIMING
struct CpamBuildTrace {
  const char* route = "unknown";
  double touch_seconds = -1.0;
  double fill_seconds = -1.0;
  double sort_seconds = 0.0;
  bool has_separate_fill = false;
  size_t input_size = 0;
};

inline double cpam_now_seconds() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

inline CpamBuildTrace& cpam_build_trace() {
  static thread_local CpamBuildTrace trace;
  return trace;
}

inline void cpam_set_build_trace(const char* route, double touch_seconds,
                                 double fill_seconds, double sort_seconds,
                                 bool has_separate_fill, size_t input_size) {
  cpam_build_trace() = CpamBuildTrace{
      .route = route,
      .touch_seconds = touch_seconds,
      .fill_seconds = fill_seconds,
      .sort_seconds = sort_seconds,
      .has_separate_fill = has_separate_fill,
      .input_size = input_size,
  };
}
#else
struct CpamBuildTrace {
  const char* route = "timing-disabled";
  double touch_seconds = -1.0;
  double fill_seconds = -1.0;
  double sort_seconds = 0.0;
  bool has_separate_fill = false;
  size_t input_size = 0;
};

inline double cpam_now_seconds() { return 0.0; }

inline CpamBuildTrace& cpam_build_trace() {
  static thread_local CpamBuildTrace trace;
  return trace;
}

inline void cpam_set_build_trace([[maybe_unused]] const char* route,
                                 [[maybe_unused]] double touch_seconds,
                                 [[maybe_unused]] double fill_seconds,
                                 [[maybe_unused]] double sort_seconds,
                                 [[maybe_unused]] bool has_separate_fill,
                                 [[maybe_unused]] size_t input_size) {}
#endif

template <typename T>
inline void cpam_parallel_first_touch(T* buffer, size_t n) {
  if (buffer == nullptr || n == 0) return;

  long ps = ::sysconf(_SC_PAGESIZE);
  const size_t page_size = (ps > 0) ? static_cast<size_t>(ps) : 4096ull;
  const size_t total_bytes = n * sizeof(T);
  const size_t pages = (total_bytes + page_size - 1) / page_size;
  const size_t workers = std::max<size_t>(1, static_cast<size_t>(parlay::num_workers()));
  auto* bytes = reinterpret_cast<unsigned char*>(buffer);

  parlay::parallel_for(
      0, workers,
      [&](size_t w) {
        size_t p0 = (pages * w) / workers;
        size_t p1 = (pages * (w + 1)) / workers;
        for (size_t p = p0; p < p1; ++p) {
          size_t offset = p * page_size;
          if (offset < total_bytes) bytes[offset] = 0;
        }
      },
      1);
}

template <typename T>
struct is_std_pair : std::false_type {};

template <typename A, typename B>
struct is_std_pair<std::pair<A, B>> : std::true_type {};

template <typename T>
inline constexpr bool is_std_pair_v = is_std_pair<T>::value;

template <typename T>
inline constexpr bool simd_sample_sort_pair_compatible_v = false;

template <typename K, typename V>
inline constexpr bool simd_sample_sort_pair_compatible_v<std::pair<K, V>> =
    std::is_pointer_v<V>;

template <typename K, typename V>
inline constexpr bool simd_sample_sort_pair_compatible_v<parlay::KVPair<K, V>> =
    std::is_integral_v<V>;

template <typename T>
concept HasCurveCodeMember = requires(T t) {
  typename T::CurveCode;
  { t.code } -> std::convertible_to<typename T::CurveCode>;
  { t.SetMember(std::declval<typename T::CurveCode const&>()) };
};

template <typename T>
struct is_materializable_curve_pair : std::false_type {};

template <typename K, typename V>
struct is_materializable_curve_pair<std::pair<K, V>>
    : std::bool_constant<std::is_pointer_v<V> && HasCurveCodeMember<K>> {};

template <typename K, typename V>
struct is_materializable_curve_pair<parlay::KVPair<K, V>>
    : std::bool_constant<std::is_integral_v<V> &&
                         std::is_same_v<K, uint64_t>> {};

template <typename T>
inline constexpr bool is_materializable_curve_pair_v =
    is_materializable_curve_pair<T>::value;

template <typename T>
inline constexpr bool dependent_false_v = false;

template <typename OutputType, typename InputType>
inline OutputType make_cpam_output_entry(InputType& in) {
  if constexpr (std::same_as<OutputType, InputType>) {
    return in;
  } else if constexpr (std::same_as<OutputType, typename InputType::AT>) {
    return in.aug;
  } else if constexpr (requires {
                         typename OutputType::first_type;
                         typename OutputType::second_type;
                       } &&
                       std::same_as<typename OutputType::first_type,
                                    typename InputType::AT>) {
    return OutputType{in.aug, &in};
  } else if constexpr (requires {
                         typename OutputType::first_type;
                         typename OutputType::second_type;
                       } &&
                       std::is_integral_v<typename OutputType::first_type> &&
                       std::is_integral_v<typename OutputType::second_type>) {
    return OutputType{
        static_cast<typename OutputType::first_type>(in.aug.code),
        static_cast<typename OutputType::second_type>(
            reinterpret_cast<uintptr_t>(&in))};
  } else {
    static_assert(dependent_false_v<OutputType>, "Unsupported output type");
  }
}

template <typename Entry, typename InputType>
inline void cpam_set_curve_entry_payload(Entry& entry, InputType* payload) {
  if constexpr (std::is_pointer_v<typename Entry::second_type>) {
    entry.second = payload;
  } else if constexpr (std::is_integral_v<typename Entry::second_type>) {
    entry.second = static_cast<typename Entry::second_type>(
        reinterpret_cast<uintptr_t>(payload));
  } else {
    static_assert(dependent_false_v<Entry>, "Unsupported curve entry payload");
  }
}

template <typename Entry, typename Code>
inline void cpam_set_curve_entry_key(Entry& entry, Code const& code) {
  if constexpr (std::is_integral_v<typename Entry::first_type>) {
    entry.first = static_cast<typename Entry::first_type>(code);
  } else if constexpr (requires { entry.first.SetMember(code); }) {
    entry.first.SetMember(code);
  } else {
    static_assert(dependent_false_v<Entry>, "Unsupported curve entry key");
  }
}

#if PSI_CPAM_HAS_SIMD_SAMPLE_SORT
template <typename filling_curve_t, typename key_entry_pointer,
          typename InputIterator>
auto cpam_sample_sort_simd_pair(slice<InputIterator, InputIterator> A,
                                bool stable = false) {
  using input_value_type =
      typename slice<InputIterator, InputIterator>::value_type;
  using key_type = typename key_entry_pointer::first_type;
  using value_ptr_type = typename key_entry_pointer::second_type;
  static_assert(std::is_integral_v<key_type>,
                "SIMD CPAM adapter expects integral key storage.");
  static_assert(std::is_integral_v<value_ptr_type>,
                "SIMD CPAM adapter expects integral payload storage.");

  sequence<key_entry_pointer> tmp =
      sequence<key_entry_pointer>::uninitialized(A.size());
  double touch_start = cpam_now_seconds();
  cpam_parallel_first_touch(tmp.begin(), tmp.size());
  double touch_seconds = cpam_now_seconds() - touch_start;
  double fill_seconds = 0.0;

  struct LazyCurveCodePrepare {
    input_value_type* input_base;
    key_entry_pointer* tmp_base;

    void prepare_entry(key_entry_pointer& entry) const {
      size_t idx = static_cast<size_t>(&entry - tmp_base);
      auto* point = input_base + idx;
      auto const code = filling_curve_t::Encode(*point);
      point->SetAugMember(code);
      cpam_set_curve_entry_key(entry, code);
      cpam_set_curve_entry_payload(entry, point);
    }

    template <typename SliceT>
    void prepare_range(SliceT s) const {
      parallel_for(
          0, s.size(),
          [&](size_t i) { prepare_entry(s[i]); },
          1024);
    }
  };

  LazyCurveCodePrepare lazy_prepare{
      .input_base = A.begin(),
      .tmp_base = tmp.begin(),
  };

  double sort_start = cpam_now_seconds();
  parlay::internal::ScratchBuffer<key_entry_pointer> scratch(tmp.size());
  auto tmp_slice = parlay::make_slice(tmp.begin(), tmp.end());
  auto scratch_slice =
      parlay::make_slice(scratch.data(), scratch.data() + tmp.size());

  if (tmp.size() < (std::numeric_limits<unsigned int>::max)()) {
    if (stable) {
      parlay::internal_simd::sample_sort_inplace_<false, true, 0, unsigned int>(
          tmp_slice, scratch_slice, true, lazy_prepare);
    } else {
      parlay::internal_simd::sample_sort_inplace_<false, false, 0,
                                                  unsigned int>(
          tmp_slice, scratch_slice, true, lazy_prepare);
    }
  } else {
    if (stable) {
      parlay::internal_simd::sample_sort_inplace_<false, true, 0, size_t>(
          tmp_slice, scratch_slice, true, lazy_prepare);
    } else {
      parlay::internal_simd::sample_sort_inplace_<false, false, 0, size_t>(
          tmp_slice, scratch_slice, true, lazy_prepare);
    }
  }
  double sort_seconds = cpam_now_seconds() - sort_start;
  cpam_set_build_trace("simd-pretouch-lazybindcode-sort", touch_seconds,
                       fill_seconds, sort_seconds, true, A.size());
  return tmp;
}
#endif

template <typename filling_curve_t, typename key_entry_pointer,
          typename InputIterator>
auto cpam_sample_sort_materialized_pair(slice<InputIterator, InputIterator> A,
                                        bool stable = false) {
  using key_type = typename key_entry_pointer::first_type;
  using value_ptr_type = typename key_entry_pointer::second_type;

  static_assert(std::is_integral_v<key_type>,
                "Materialized CPAM sort expects integral key storage.");
  static_assert(std::is_integral_v<value_ptr_type>,
                "Materialized CPAM sort expects integral payload storage.");

  double fill_start = cpam_now_seconds();
  sequence<key_entry_pointer> tmp =
      sequence<key_entry_pointer>::uninitialized(A.size());
  parallel_for(0, A.size(), [&](size_t i) {
    auto const code = filling_curve_t::Encode(A[i]);
    A[i].SetAugMember(code);
    tmp[i] = key_entry_pointer{
        static_cast<key_type>(code),
        static_cast<value_ptr_type>(reinterpret_cast<uintptr_t>(&A[i]))};
  });
  double fill_seconds = cpam_now_seconds() - fill_start;

  auto less_by_code = [&](auto const& a, auto const& b) {
    return a.first < b.first;
  };
  auto tmp_slice = parlay::make_slice(tmp.begin(), tmp.end());
  double sort_start = cpam_now_seconds();
  if (stable) {
    bucket_sort(tmp_slice, less_by_code, true);
  } else {
    quicksort(tmp.begin(), tmp.size(), less_by_code);
  }
  double sort_seconds = cpam_now_seconds() - sort_start;
  cpam_set_build_trace("scalar-precompute-fill-sort", -1.0, fill_seconds,
                       sort_seconds, true, A.size());
  return tmp;
}

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
    Out[j] = make_cpam_output_entry<OutputType>(In[j]);
  }
  // for (size_t j = 0; j < l; j++) {
  //   std::cout << In[j].aug.code << " ";
  // }

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
      return make_cpam_output_entry<output_value_type>(*pt);
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
  (void)sizeof(input_value_type);
  (void)less;

#if PSI_CPAM_HAS_SIMD_SAMPLE_SORT
  if constexpr (is_materializable_curve_pair_v<key_entry_pointer> &&
                simd_sample_sort_pair_compatible_v<key_entry_pointer>) {
    return cpam_sample_sort_simd_pair<filling_curve_t, key_entry_pointer>(
        A, stable);
  }
#endif

#if PSI_CPAM_PRECOMPUTE_CODE_SORT
  if constexpr (is_materializable_curve_pair_v<key_entry_pointer>) {
    return cpam_sample_sort_materialized_pair<filling_curve_t,
                                              key_entry_pointer>(A, stable);
  }
#endif

  sequence<key_entry_pointer> R =
      sequence<key_entry_pointer>::uninitialized(A.size());
  double sort_start = cpam_now_seconds();
  if (A.size() < (std::numeric_limits<unsigned int>::max)()) {
    sample_sort_<filling_curve_t, unsigned int>(A, make_slice(R), less, stable);
  } else {
    sample_sort_<filling_curve_t, size_t>(A, make_slice(R), less, stable);
  }
  double sort_seconds = cpam_now_seconds() - sort_start;
  cpam_set_build_trace("legacy-interleaved-fill-sort", -1.0, -1.0,
                       sort_seconds, false, A.size());
  return R;
}

}  // namespace cpam
}  // namespace internal

}  // namespace parlay
#endif  // CPAM_SAMPLE_SORT_H
