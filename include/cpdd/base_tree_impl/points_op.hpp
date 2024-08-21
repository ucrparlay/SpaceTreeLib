#pragma once

#include "../base_tree.h"
#include <algorithm>
#include <utility>

namespace cpdd {
template<typename Point, uint_fast8_t kBDO>
inline void BaseTree<Point, kBDO>::SamplePoints(Slice In, Points& arr) {
    auto size = arr.size();
    auto n = In.size();
    auto indexs = parlay::sequence<uint64_t>::uninitialized(size);
    for (size_t i = 0; i < size; i++) {
        indexs[i] = parlay::hash64(i) % n;
    }
    std::sort(indexs.begin(), indexs.end());
    for (size_t i = 0; i < size; i++) {
        arr[i] = In[indexs[i]];
    }
    return;
}

template<typename Point, uint_fast8_t kBDO>
inline typename BaseTree<Point, kBDO>::BucketType
BaseTree<Point, kBDO>::FindBucket(const Point& p, const HyperPlaneSeq& pivots) {
    BucketType k(1);
    while (k <= kPivotNum) {
        k = k * 2 + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[pivots[k].second], pivots[k].first));
    }
    assert(pivots[k].first == -1);
    return pivots[k].second;
}

template<typename Point, uint_fast8_t kBDO>
void BaseTree<Point, kBDO>::Partition(Slice A, Slice B, const size_t n,
                                      const HyperPlaneSeq& pivots,
                                      parlay::sequence<BallsType>& sums) {
    size_t num_block = (n + kBlockSize - 1) >> kLog2Base;
    parlay::sequence<parlay::sequence<BallsType>> offset(
        num_block, parlay::sequence<BallsType>(kBucketNum));
    assert(offset.size() == num_block && offset[0].size() == kBucketNum &&
           offset[0][0] == 0);
    parlay::parallel_for(0, num_block, [&](size_t i) {
        for (size_t j = i << kLog2Base; j < std::min((i + 1) << kLog2Base, n);
             j++) {
            DimsType k = FindBucket(A[j], pivots);
            offset[i][k]++;
            /*offset[i][std::move(find_bucket(A[j], pivots))]++;*/
        }
    });

    sums = parlay::sequence<BallsType>(kBucketNum);
    for (size_t i = 0; i < num_block; i++) {
        auto t = offset[i];
        offset[i] = sums;
        for (BucketType j = 0; j < kBucketNum; ++j) {
            sums[j] += t[j];
        }
    }

    parlay::parallel_for(0, num_block, [&](size_t i) {
        auto v = parlay::sequence<BallsType>::uninitialized(kBucketNum);
        size_t tot = 0, s_offset = 0;
        for (BucketType k = 0; k < kBucketNum - 1; ++k) {
            v[k] = tot + offset[i][k];
            tot += sums[k];
            s_offset += offset[i][k];
        }
        v[kBucketNum - 1] = tot + ((i << kLog2Base) - s_offset);
        for (size_t j = i << kLog2Base; j < std::min((i + 1) << kLog2Base, n);
             j++) {
            DimsType k = FindBucket(A[j], pivots);
            // B[v[std::move( find_bucket( A[j], pivots ) )]++] = A[j];
            B[v[k]++] = A[j];
        }
    });

    return;
}

template<typename Point, uint_fast8_t kBDO>
typename BaseTree<Point, kBDO>::PointsIter
BaseTree<Point, kBDO>::SerialPartition(Slice In, DimsType d) {
    size_t n = In.size();
    std::ranges::nth_element(In.begin(), In.begin() + n / 2, In.end(),
                             [&](const Point& p1, const Point& p2) {
                                 return Num::Lt(p1.pnt[d], p2.pnt[d]);
                             });

    std::ranges::subrange _2ndGroup = std::ranges::partition(
        In.begin(), In.begin() + n / 2,
        [&](const Point& p) { return Num::Lt(p.pnt[d], In[n / 2].pnt[d]); });

    if (_2ndGroup.begin() == In.begin()) {  // NOTE: handle duplicated medians
        _2ndGroup = std::ranges::partition(
            In.begin() + n / 2, In.end(), [&](const Point& p) {
                return Num::Eq(p.pnt[d], In[n / 2].pnt[d]);
            });  // NOTE: now all duplicated median is on the left
    }
    return _2ndGroup.begin();
}

template<typename Point, uint_fast8_t kBDO>
template<typename Leaf, typename Interior, bool granularity>
void BaseTree<Point, kBDO>::PrepareRebuild(Node* T, Slice In, Points& wx,
                                           Points& wo) {
    // TODO: add dispatch tag
    wo = Points::uninitialized(T->size + In.size());
    wx = Points::uninitialized(T->size + In.size());
    parlay::parallel_for(0, In.size(), [&](size_t j) { wx[j] = In[j]; });
    FlattenRec<Leaf, Interior, Slice, granularity>(
        T, wx.cut(In.size(), wx.size()));
    DeleteTreeRecursive<Leaf, Interior, granularity>(T);
    return;
}

// NOTE: retrive the bucket tag of Point p from the skeleton tags
template<typename Point, uint_fast8_t kBDO>
template<IsBinaryNode Interior>
typename BaseTree<Point, kBDO>::BucketType BaseTree<Point, kBDO>::RetriveTag(
    const Point& p, const NodeTagSeq& tags) {
    BucketType k(1);
    while (k <= kPivotNum && (!tags[k].first->is_leaf)) {
        k = k * 2 + 1 -
            static_cast<BucketType>(Num::Lt(
                p.pnt[static_cast<Interior*>(tags[k].first)->split.second],
                static_cast<Interior*>(tags[k].first)->split.first));
    }
    assert(tags[k].second < kBucketNum);
    return tags[k].second;
}

template<typename Point, uint_fast8_t kBDO>
template<IsMultiNode Interior>
typename BaseTree<Point, kBDO>::BucketType BaseTree<Point, kBDO>::RetriveTag(
    const Point& p, const NodeTagSeq& tags) {
    BucketType k(1);
    while (k <= kPivotNum && (!tags[k].first->is_leaf)) {
        k = static_cast<Interior*>(tags[k].first)->SeievePoint(p, k);
    }
    assert(tags[k].second < kBucketNum);
    return tags[k].second;
}

// NOTE: seieve Points from range A to range B, using the skeleton tags. The
// sums is the number of elemenets within each bucket, the tags_num is the total
// number of buckets in the skeleton
template<typename Point, uint_fast8_t kBDO>
template<typename Interior>
void BaseTree<Point, kBDO>::SeievePoints(Slice A, Slice B, const size_t n,
                                         const NodeTagSeq& tags,
                                         parlay::sequence<BallsType>& sums,
                                         const BucketType tags_num) {
    size_t num_block = (n + kBlockSize - 1) >> kLog2Base;
    parlay::sequence<parlay::sequence<BallsType>> offset(
        num_block, parlay::sequence<BallsType>(tags_num));
    assert(offset.size() == num_block && offset[0].size() == tags_num &&
           offset[0][0] == 0);
    parlay::parallel_for(0, num_block, [&](size_t i) {
        for (size_t j = i << kLog2Base; j < std::min((i + 1) << kLog2Base, n);
             j++) {
            auto k = RetriveTag<Interior>(A[j], tags);
            offset[i][k]++;
        }
    });

    sums = parlay::sequence<BallsType>(tags_num);
    for (size_t i = 0; i < num_block; i++) {
        auto t = offset[i];
        offset[i] = sums;
        for (int j = 0; j < tags_num; j++) {
            sums[j] += t[j];
        }
    }

    parlay::parallel_for(0, num_block, [&](size_t i) {
        auto v = parlay::sequence<BallsType>::uninitialized(tags_num);
        int tot = 0, s_offset = 0;
        for (int k = 0; k < tags_num - 1; k++) {
            v[k] = tot + offset[i][k];
            tot += sums[k];
            s_offset += offset[i][k];
        }
        v[tags_num - 1] = tot + ((i << kLog2Base) - s_offset);
        for (size_t j = i << kLog2Base; j < std::min((i + 1) << kLog2Base, n);
             j++) {
            auto k = RetriveTag<Interior>(A[j], tags);
            B[v[k]++] = A[j];
        }
    });

    return;
}
}  // namespace cpdd
