#pragma once

#include "../base_tree.h"
#include <algorithm>
#include <utility>

namespace cpdd {
template<typename Point, uint8_t kBDO>
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

template<typename Point, uint8_t kBDO>
inline uint_fast8_t BaseTree<Point, kBDO>::FindBucket(
    const Point& p, const HyperPlaneSeq& pivots) {
    uint_fast8_t k = 1;
    while (k <= kPivotNum) {
        // TODO: remove conditional judge
        k = Num::Lt(p.pnt[pivots[k].second], pivots[k].first) ? k << 1
                                                              : k << 1 | 1;
    }
    assert(pivots[k].first == -1);
    return pivots[k].second;
}

template<typename Point, uint8_t kBDO>
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

template<typename Point, uint8_t kBDO>
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

}  // namespace cpdd
