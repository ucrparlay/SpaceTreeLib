#pragma once

#include <cstddef>
#include <cstdlib>

// NOTE: Coordinate and typename definitions
#ifdef CCP
using Coord = long;
// using Coord = double;
#else
using Coord = unsigned long long;
// using Coord = double;
#endif  // CCP

using Typename = Coord;

// NOTE: KNN size
static constexpr double kBatchQueryRatio = 0.01;
static constexpr size_t kBatchQueryOsmSize = 10000000;

// NOTE: rectangle numbers
static constexpr int kRangeQueryNum = 50000;
static constexpr int kSmallRangeQueryNum = 100000000;
static constexpr int kMediumRangeQueryNum = 1000000;
static constexpr int kLargeRangeQueryNum = 100000;
static constexpr int kSingleQueryLogRepeatNum = 100;

// NOTE: rectangle numbers for inba ratio
static constexpr int rangeQueryNumInbaRatio = 50000;

// NOTE: insert batch ratio for inba ratio
static constexpr double insertBatchInbaRatio = 0.001;

// NOTE: knn batch ratio for inba ratio
static constexpr double knnBatchInbaRatio = 0.1;

// NOTE: Insert Ratio when summary
static constexpr double kBatchInsertRatio = 0.01;

// NOTE: Diff Ratio when summary
static constexpr double kBatchDiffTotalRatio = 0.01;
static constexpr double kBatchDiffOverlapRatio = 0.2;

// NOTE: rectangle type used in summary
static constexpr int summaryRangeQueryType = 2;

// NOTE: range query num in summary
static constexpr int kSummaryRangeQueryNum = 50000;

// NOTE: helper for delete type
enum DeleteType { kBatchDelete, kBatchDiff };

// * [a,b)
inline size_t get_random_index(size_t a, size_t b, [[maybe_unused]] int seed) {
  return size_t((rand() % (b - a)) + a);
}
