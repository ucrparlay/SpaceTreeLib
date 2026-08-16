#ifndef PSI_BASE_TREE_IMPL_POINTS_OP_HPP_
#define PSI_BASE_TREE_IMPL_POINTS_OP_HPP_

#include <algorithm>
#include <utility>

#include "../base_tree.h"

namespace psi
{
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::sample_points(
	slice_type in, points_type &arr)
{
	auto size = arr.size();
	auto n = in.size();
	auto indexs = parlay::sequence<uint64_t>::uninitialized(size);
	for (size_t i = 0; i < size; i++) {
		indexs[i] = parlay::hash64(i) % n;
	}
	std::ranges::sort(indexs);
	for (size_t i = 0; i < size; i++) {
		arr[i] = in[indexs[i]];
	}
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::bucket_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::find_bucket(
	Point const &p, hyper_plane_seq_type const &pivots)
{
	bucket_type k(1);
	while (k <= pivot_num) {
		k = k * 2 + 1 -
		    static_cast<bucket_type>(num_type::lt(
			    p.pnt[pivots[k].second], pivots[k].first));
	}
	assert(pivots[k].first == 0);
	return pivots[k].second;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::partition(
	slice_type A, slice_type B, size_t const n,
	hyper_plane_seq_type const &pivots, parlay::sequence<balls_type> &sums)
{
	size_t num_block = (n + block_size - 1) >> log2_base;
	parlay::sequence<parlay::sequence<balls_type>> offset(
		num_block, parlay::sequence<balls_type>(bucket_num));
	assert(offset.size() == num_block);
	assert(offset[0].size() == bucket_num);
	assert(offset[0][0] == 0);
	parlay::parallel_for(0, num_block, [&](size_t i) {
		for (size_t j = i << log2_base;
		     j < std::min((i + 1) << log2_base, n); j++) {
			dims_type k = find_bucket(A[j], pivots);
			offset[i][k]++;
			/*offset[i][std::move(find_bucket(A[j], pivots))]++;*/
		}
	});

	sums = parlay::sequence<balls_type>(bucket_num);
	for (size_t i = 0; i < num_block; i++) {
		auto t = offset[i];
		offset[i] = sums;
		for (bucket_type j = 0; j < bucket_num; ++j) {
			sums[j] += t[j];
		}
	}

	parlay::parallel_for(0, num_block, [&](size_t i) {
		auto v =
			parlay::sequence<balls_type>::uninitialized(bucket_num);
		size_t tot = 0, s_offset = 0;
		for (bucket_type k = 0; k < bucket_num - 1; ++k) {
			v[k] = tot + offset[i][k];
			tot += sums[k];
			s_offset += offset[i][k];
		}
		v[bucket_num - 1] = tot + ((i << log2_base) - s_offset);
		for (size_t j = i << log2_base;
		     j < std::min((i + 1) << log2_base, n); j++) {
			dims_type k = find_bucket(A[j], pivots);
			// B[v[std::move( find_bucket( A[j], pivots ) )]++] =
			// A[j];
			B[v[k]++] = A[j];
		}
	});

	return;
}

// NOTE: retrieve the bucket tag of Point p from the skeleton tags
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <is_binary_node interior_type>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::bucket_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::retrieve_tag(
	Point const &p, node_tag_seq_type const &tags)
{
	bucket_type k(1);
	while (k <= pivot_num && (!tags[k].first->is_leaf)) {
		k = k * 2 + 1 -
		    static_cast<bucket_type>(num_type::lt(
			    p.pnt[static_cast<interior_type *>(tags[k].first)
					  ->split.second],
			    static_cast<interior_type *>(tags[k].first)
				    ->split.first));
	}
	assert(tags[k].second < bucket_num);
	return tags[k].second;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <is_multi_node interior_type>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::bucket_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::retrieve_tag(
	Point const &p, node_tag_seq_type const &tags)
{
	bucket_type k(1);
	while (k <= pivot_num && (!tags[k].first->is_leaf)) {
		k = static_cast<interior_type *>(tags[k].first)
			    ->sieve_point(p, k);
	}
	assert(tags[k].second < bucket_num);
	return tags[k].second;
}

// NOTE: sieve points_type from range A to range B, using the skeleton tags. The
// sums is the number of elemenets within each bucket, the tags_num is the total
// number of buckets in the skeleton
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename interior_type>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::sieve_points(
	slice_type A, slice_type B, size_t const n,
	node_tag_seq_type const &tags, parlay::sequence<balls_type> &sums,
	bucket_type const tags_num)
{
	size_t num_block = (n + block_size - 1) >> log2_base;
	parlay::sequence<parlay::sequence<balls_type>> offset(
		num_block, parlay::sequence<balls_type>(tags_num));
	assert(offset.size() == num_block && offset[0].size() == tags_num &&
	       offset[0][0] == 0);
	parlay::parallel_for(0, num_block, [&](size_t i) {
		for (size_t j = i << log2_base;
		     j < std::min((i + 1) << log2_base, n); j++) {
			auto k = retrieve_tag<interior_type>(A[j], tags);
			offset[i][k]++;
		}
	});

	sums = parlay::sequence<balls_type>(tags_num);
	for (size_t i = 0; i < num_block; i++) {
		auto t = offset[i];
		offset[i] = sums;
		for (bucket_type j = 0; j < tags_num; ++j) {
			sums[j] += t[j];
		}
	}

	parlay::parallel_for(0, num_block, [&](size_t i) {
		auto v = parlay::sequence<balls_type>::uninitialized(tags_num);
		size_t tot = 0, s_offset = 0;
		for (bucket_type k = 0; k < tags_num - 1; ++k) {
			v[k] = tot + offset[i][k];
			tot += sums[k];
			s_offset += offset[i][k];
		}
		v[tags_num - 1] = tot + ((i << log2_base) - s_offset);
		for (size_t j = i << log2_base;
		     j < std::min((i + 1) << log2_base, n); j++) {
			auto k = retrieve_tag<interior_type>(A[j], tags);
			B[v[k]++] = A[j];
		}
	});

	return;
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_POINTS_OP_HPP_
