
#ifndef DOVETAIL_INTERNAL_INTEGER_SORT_H_
#define DOVETAIL_INTERNAL_INTEGER_SORT_H_

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "dovetail_counting_sort.h"
#include "parlay/delayed_sequence.h"
#include "parlay/internal/get_time.h"
#include "parlay/internal/sequence_ops.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/monoid.h"
#include "parlay/parallel.h"
#include "parlay/range.h"
#include "parlay/relocation.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/utilities.h"

namespace parlay
{
namespace internal
{
namespace cpam
{

constexpr size_t radix = 8;
constexpr size_t max_buckets = 1 << radix;

// Use a smaller base case threshold for debugging so that test
// cases do not need to use extremely large sequences in order
// to achieve adequate coverage.
#ifdef DEBUG
constexpr size_t parlay_integer_sort_base_case_size = 128;
#else
constexpr size_t parlay_integer_sort_base_case_size = 1 << 17;
#endif

// a bottom up radix sort
template <typename InIterator, typename OutIterator, class GetKey>
void seq_radix_sort_(slice<InIterator, InIterator> in,
		     slice<OutIterator, OutIterator> out, GetKey const &g,
		     size_t bits, bool inplace = true)
{
	size_t n = in.size();
	if (n == 0)
		return;
	size_t counts[max_buckets + 1];
	bool swapped = false;
	size_t bit_offset = 0;
	while (bits > 0) {
		size_t round_bits = (std::min)(radix, bits);
		size_t num_buckets = (size_t{1} << round_bits);
		size_t mask = num_buckets - 1;

		if (swapped) {
			auto get_key = [&](size_t i) -> size_t {
				return (g(out[i]) >> bit_offset) & mask;
			};
			seq_count_sort_<uninitialized_relocate_tag>(
				out, in, delayed_seq<size_t>(n, get_key),
				counts, num_buckets);
		}

		else {
			auto get_key = [&](size_t i) -> size_t {
				return (g(in[i]) >> bit_offset) & mask;
			};

			seq_count_sort_<uninitialized_relocate_tag>(
				in, out, delayed_seq<size_t>(n, get_key),
				counts, num_buckets);
		}

		bits = bits - round_bits;
		bit_offset = bit_offset + round_bits;
		swapped = !swapped;
	}

	if (swapped && inplace) {
		parlay::uninitialized_relocate(out.begin(), out.end(),
					       in.begin());
	} else if (!swapped && !inplace) {
		parlay::uninitialized_relocate(in.begin(), in.end(),
					       out.begin());
	}
}

// wrapper to reduce copies and avoid modifying in when not inplace
// in and tmp can be the same, but out must be different
template <typename inplace_tag, typename assignment_tag, typename InIterator,
	  typename OutIterator, typename TmpIterator, class GetKey>
void seq_radix_sort(slice<InIterator, InIterator> in,
		    slice<OutIterator, OutIterator> out,
		    [[maybe_unused]] slice<TmpIterator, TmpIterator> tmp,
		    GetKey const &g, size_t key_bits)
{
	if constexpr (inplace_tag::value == true) {
		static_assert(std::is_same_v<assignment_tag,
					     uninitialized_relocate_tag>);
		seq_radix_sort_(in, out, g, key_bits, true);
	} else {
		bool odd = ((key_bits - 1) / radix) & 1;
		size_t n = in.size();
		if (odd) {
			// We could just use assign_dispatch(tmp[i], in[i]) for
			// each i, but we can optimize better by calling
			// uninitialized_relocate, since this has the ability to
			// memcpy multiple elements at once
			if constexpr (std::is_same_v<
					      assignment_tag,
					      uninitialized_relocate_tag>) {
				parlay::uninitialized_relocate(
					in.begin(), in.end(), tmp.begin());
			} else {
				for (size_t i = 0; i < n; i++)
					assign_uninitialized(tmp[i], in[i]);
			}
			seq_radix_sort_(tmp, out, g, key_bits, false);
		} else {
			if constexpr (std::is_same_v<
					      assignment_tag,
					      uninitialized_relocate_tag>) {
				parlay::uninitialized_relocate(
					in.begin(), in.end(), out.begin());
			} else {
				for (size_t i = 0; i < n; i++)
					assign_uninitialized(out[i], in[i]);
			}
			seq_radix_sort_(out, tmp, g, key_bits, true);
		}
	}
}

// a top down recursive radix sort
// g extracts the integer keys from in
// key_bits specifies how many bits there are left
//
// inplace_tag must be one of std::true_type or std::false_type. If inplace_tag
// is true_type, then the output of the sort will remain in in. If inplace_tag
// is false, then the output of the sort will be written into out. Furthermore,
// if inplace_tag is std::true_type, then tmp must point to the same range as in
//
// assignment_tag must be one of uninitialized_copy_tag or
// uninitialized_relocate_tag. This indicates the manner in which the input is
// moved from in to out. If inplace_tag is std::true_type, then this must always
// be uninitialized_relocate_tag. If inplace_tag is std::false_type, then
// assignment_tag can be either uninitialized_copy_tag, if the input is to be
// copied from in to out, or uninitialized_relocate_tag if the input is to be
// (destructively) moved from in to out.
//
template <typename inplace_tag, typename assignment_tag, typename InIterator,
	  typename OutIterator, typename TmpIterator, typename Get_Key>
sequence<size_t> integer_sort_r(slice<InIterator, InIterator> in,
				slice<OutIterator, OutIterator> out,
				slice<TmpIterator, TmpIterator> tmp,
				Get_Key const &g, size_t key_bits,
				size_t num_buckets, float parallelism = 1.0)
{
	// Parameter in
	//  [Preconditions]
	//  - in owns valid initialized memory
	//  - If inplace_tag is std::true_type, then in points to mutable
	//  (non-const) memory
	//  - If assignment_tag is uninitialized_copy_tag, then in points to
	//  objects of a type
	//    that is copy constructible
	//  [Postconditions]
	//  - If inplace_tag is std::true_type, then in contains the sorted
	//  results
	//  - If inplace_tag is std::false_type and assignment_tag is
	//  uninitialized_copy_tag,
	//    then in contains the original unmodified input
	//  - If inplace_tag is std::false_type and assignment_tag is
	//  uninitialized_relocate_tag, then
	//    in points to uninitialized memory
	//
	// Parameter out
	//  [Preconditions]
	//  - out points to uninitialized memory
	//  [Postconditions]
	//  - If inplace_tag is std::true_type, then out points to uninitialized
	//  memory
	//  - If inplace_tag is std::false_type, then out points to the
	//  resulting sorted objects
	//
	// Parameter tmp
	//  [Preconditions]
	//  - If inplace_tag is std::true_type, then tmp points to the same
	//  range as in
	//  - If inplace_tag is std::false_type and assignment_tag is
	//  uninitialized_copy_tag,
	//    then tmp points to uninitialized memory
	//  - If inplace_tag is std::false_type and assignment_tag is
	//  uninitialized_relocate_tag,
	//    then tmp points to the same range as in
	//  [Postconditions]
	//  - If tmp does not point to the same range as in, then it points to
	//  uninitialized memory

	// Ranges should all have the same underlying type
	static_assert(std::is_same_v<
		      range_value_type_t<slice<InIterator, InIterator>>,
		      range_value_type_t<slice<OutIterator, OutIterator>>>);
	static_assert(std::is_same_v<
		      range_value_type_t<slice<InIterator, InIterator>>,
		      range_value_type_t<slice<TmpIterator, TmpIterator>>>);

	// assignment_type can only be one of uninitialized_copy_tag or
	// uninitialized_relocate_tag
	static_assert(
		std::is_same_v<assignment_tag, uninitialized_copy_tag> ||
		std::is_same_v<assignment_tag, uninitialized_relocate_tag>);

	// assignment_tag is only allowed to be uninitialized_copy_tag if
	// inplace_tag is std::false type
	static_assert(
		inplace_tag::value == false ||
		std::is_same_v<assignment_tag, uninitialized_relocate_tag>);

	using T = typename slice<InIterator, InIterator>::value_type;
	timer t("integer sort", false);

	size_t n = in.size();
	size_t cache_per_thread = 1000000;
	auto sz = 2 * (size_t)sizeof(T) * n / cache_per_thread;
	size_t base_bits = sz > 0 ? log2_up(sz) : 0;
	// keep between 8 and 13
	base_bits = std::max<size_t>(8, std::min<size_t>(13, base_bits));
	sequence<size_t> offsets;
	bool one_bucket;
	bool return_offsets = (num_buckets > 0);

	if (key_bits == 0) {
		// If the sort is not inplace, the final result needs to
		// be moved into out since it is currently in in.
		if constexpr (inplace_tag::value == false) {
			// We could just do a parallel for loop and copy/move
			// the elements from in to out using
			// assign_dispatch(out[i], in[i], assignment_tag()), but
			// this would not allow us to take advantage of the
			// possibly more efficient uninitialized_relocate_n,
			// which can memcpy multiple elements at a time to save
			// on performing every copy individually.
			if constexpr (std::is_same_v<
					      assignment_tag,
					      uninitialized_relocate_tag>) {
				parlay::uninitialized_relocate(
					in.begin(), in.end(), out.begin());
			} else {
				parallel_for(0, in.size(), [&](size_t i) {
					assign_uninitialized(
						out[i],
						in[i]); // Copy from in[i] to
							// out[i]
				});
			}
		}
		return sequence<size_t>();
	}
	// for small inputs or little parallelism use sequential radix sort
	else if ((n < parlay_integer_sort_base_case_size ||
		  parallelism < .0001) &&
		 !return_offsets) {
		seq_radix_sort<inplace_tag, assignment_tag>(in, out, tmp, g,
							    key_bits);
		return sequence<size_t>();
	}
	// few bits, just do a single parallel count sort
	else if (key_bits <= base_bits) {
		size_t mask = (1 << key_bits) - 1;
		auto f = [&](size_t i) {
			return static_cast<size_t>(g(in[i]) & mask);
		};
		auto get_bits = delayed_seq<size_t>(n, f);
		size_t num_bkts = (num_buckets == 0) ? (size_t{1} << key_bits)
						     : num_buckets;

		// only uses one bucket optimization (last argument) if inplace
		std::tie(offsets, one_bucket) = count_sort<assignment_tag>(
			in, out, make_slice(get_bits), num_bkts, parallelism,
			inplace_tag::value);
		t.next("count sort");

		if constexpr (inplace_tag::value == true) {
			if (!one_bucket) {
				parlay::uninitialized_relocate(
					out.begin(), out.end(), in.begin());
			}
		}

		if (return_offsets)
			return offsets;
		else
			return sequence<size_t>();

		// recursive case
	} else {
		size_t bits = 8;
		size_t shift_bits = key_bits - bits;
		size_t num_outer_buckets = (size_t{1} << bits);
		size_t num_inner_buckets =
			return_offsets ? ((size_t)1 << shift_bits) : 0;
		size_t mask = num_outer_buckets - 1;
		auto f = [&](size_t i) {
			return static_cast<size_t>((g(in[i]) >> shift_bits) &
						   mask);
		};
		auto get_bits = delayed_seq<size_t>(n, f);

		// divide into buckets
		std::tie(offsets, one_bucket) = count_sort<assignment_tag>(
			in, out, make_slice(get_bits), num_outer_buckets,
			parallelism, !return_offsets);
		if (parallelism == 1.0)
			t.next("recursive count sort");

		// if all but one bucket are empty, try again on lower bits
		if (one_bucket) {
			return integer_sort_r<inplace_tag, assignment_tag>(
				in, out, tmp, g, shift_bits, 0, parallelism);
		}

		// After this point, out is guaranteed to be initialized

		sequence<size_t> inner_offsets(return_offsets ? num_buckets + 1
							      : 0);
		if (return_offsets)
			inner_offsets[num_buckets] = n;

		// recursively sort each bucket
		parallel_for(
			0, num_outer_buckets,
			[&](size_t i) {
				size_t start = offsets[i];
				size_t end = offsets[i + 1];
				auto a = out.cut(start, end);
				auto b = tmp.cut(start, end);
				sequence<size_t> r;

				auto new_parallelism =
					(parallelism *
					 static_cast<float>(end - start)) /
					static_cast<float>(n + 1);
				r = integer_sort_r<typename std::negation<
							   inplace_tag>::type,
						   uninitialized_relocate_tag>(
					a, b, a, g, shift_bits,
					num_inner_buckets, new_parallelism);

				if (return_offsets) {
					size_t bstart =
						(std::min)(i * num_inner_buckets,
							   num_buckets);
					size_t bend =
						(std::min)((i +
							    1) * num_inner_buckets,
							   num_buckets);
					size_t m = (bend - bstart);
					for (size_t j = 0; j < m; j++)
						inner_offsets[bstart + j] =
							offsets[i] + r[j];
				}
			},
			1);
		return inner_offsets;
	}
}

// a top down recursive radix sort
// g extracts the integer keys from in
// if inplace is false then result will be placed in out,
//    otherwise they are placed in tmp
//    tmp and in can be the same (i.e. to do inplace set them equal)
// in is not directly modified, but can be indirectly if equal to tmp
// bits specifies how many bits there are in the key
//    if set to 0, then a max is taken over the keys to determine
// If num_buckets is non-zero then the output sequence will contain
// the offsets of each bucket (num_bucket of them)
// num_bucket must be less than or equal to 2^bits
template <typename inplace_tag, typename assignment_tag, typename InIterator,
	  typename OutIterator, typename TmpIterator, typename Get_Key>
sequence<size_t> integer_sort_(slice<InIterator, InIterator> in,
			       slice<OutIterator, OutIterator> out,
			       slice<TmpIterator, TmpIterator> tmp,
			       Get_Key const &g, size_t bits,
			       size_t num_buckets)
{
	if (bits == 0) {
		auto get_key = [&](size_t i) {
			return static_cast<size_t>(g(in[i]));
		};
		auto keys = delayed_seq<size_t>(in.size(), get_key);
		using key_type = std::invoke_result_t<
			Get_Key,
			typename slice<InIterator, InIterator>::value_type>;
		auto max_key =
			internal::reduce(make_slice(keys), maximum<key_type>());
		if (max_key == std::numeric_limits<key_type>::max()) {
			bits = sizeof(key_type) * 8;
		} else {
			bits = log2_up(max_key + 1);
		}
	}
	return integer_sort_r<inplace_tag, assignment_tag>(in, out, tmp, g,
							   bits, num_buckets);
}

template <typename Iterator, typename Get_Key>
void integer_sort_inplace(slice<Iterator, Iterator> in, Get_Key const &g,
			  size_t bits = 0)
{
	using value_type = typename slice<Iterator, Iterator>::value_type;
	auto tmp = internal::uninitialized_sequence<value_type>(in.size());
	integer_sort_<std::true_type, uninitialized_relocate_tag>(
		in, make_slice(tmp), in, g, bits, 0);
}

template <typename Iterator, typename Get_Key>
auto integer_sort(slice<Iterator, Iterator> in, Get_Key const &g,
		  size_t bits = 0)
{
	using value_type = typename slice<Iterator, Iterator>::value_type;
	auto out = sequence<value_type>::uninitialized(in.size());
	auto tmp = uninitialized_sequence<value_type>(in.size());
	integer_sort_<std::false_type, uninitialized_copy_tag>(
		in, make_slice(out), make_slice(tmp), g, bits, 0);
	return out;
}

// Given a sorted sequence of integers in the range [0,..,num_buckets)
// returns a sequence of length num_buckets+1 with the offset for the
// start of each integer.   If an integer does not appear, its offset
// will be the same as the next (i.e. offset[i+1]-offset[i] specifies
// how many i there are.
// The last element contains the size of the input.
template <typename Tint = size_t, typename Iterator, typename Get_Key>
sequence<Tint> get_counts(slice<Iterator, Iterator> in, Get_Key const &g,
			  size_t num_buckets)
{
	size_t n = in.size();
	if (n == 0) {
		return {};
	}
	sequence<Tint> starts(num_buckets, (Tint)0);
	sequence<Tint> ends(num_buckets, (Tint)0);
	parallel_for(0, n - 1, [&](size_t i) {
		if (g(in[i]) != g(in[i + 1])) {
			starts[g(in[i + 1])] = i + 1;
			ends[g(in[i])] = i + 1;
		};
	});
	ends[g(in[n - 1])] = n;
	return sequence<Tint>::from_function(
		num_buckets, [&](size_t i) { return ends[i] - starts[i]; });
}

template <typename Tint = size_t, typename Iterator, typename Get_Key>
auto integer_sort_with_counts(slice<Iterator, Iterator> in, Get_Key const &g,
			      size_t num_buckets)
{
	using T = typename slice<Iterator, Iterator>::value_type;
	if (in.size() == 0) {
		return std::make_pair(parlay::sequence<T>{},
				      parlay::sequence<Tint>(num_buckets));
	}
	assert(num_buckets > 0);
	size_t bits = log2_up(num_buckets);
	auto R = integer_sort(in, g, bits);
	return std::make_pair(std::move(R),
			      get_counts<Tint>(make_slice(R), g, num_buckets));
}

} // namespace cpam
} // namespace internal
} // namespace parlay

#endif // DOVETAIL_INTERNAL_INTEGER_SORT_H_
