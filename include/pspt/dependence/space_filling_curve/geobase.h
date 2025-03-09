#pragma once
#include <bits/stdc++.h>
#include <parlay/internal/binary_search.h>

#include "hilbert.h"

namespace geobase {
using namespace std;
using FT = double;
// using FT = float;
constexpr FT FT_INF_MIN = numeric_limits<FT>::min();
constexpr FT FT_INF_MAX = numeric_limits<FT>::max();
constexpr FT FT_EPS = numeric_limits<FT>::epsilon();

struct break_down {
  FT sort_time = 0;
  FT leaf_time = 0;
  FT inte_time = 0;
  FT split_time = 0;
  FT slice_time = 0;
  void clear() {
    sort_time = 0;
    leaf_time = 0;
    inte_time = 0;
    split_time = 0;
    slice_time = 0;
  }
};

inline int dcmp(const FT& x) {
  if (fabs(x) < FT_EPS) return 0;
  return x < 0 ? -1 : 1;
}

bool less_msb(unsigned int x, unsigned int y) { return x < y && x < (x ^ y); }

struct Point {
  size_t id;
  FT x, y;
  unsigned long long morton_id;
  Point() {}
  Point(FT _x, FT _y) : x(_x), y(_y) {}
  Point(unsigned _id, FT _x, FT _y) : id(_id), x(_x), y(_y) {
    // morton_id = mortonIndex();
    // morton_id = interleave_bits();
  }

  bool operator==(Point const& p) const {
    return !(id - p.id) && !dcmp(x - p.x) && !dcmp(y - p.y);
  }

  friend std::ostream& operator<<(std::ostream& os, Point const& p) {
    os << p.id << ": (" << p.x << ", " << p.y << ")";
    // os << "(" << p.x << ", " << p.y << ")";
    return os;
  }

  // bool operator < (const Point &b) const{
  //     return (morton_id < b.morton_id) ||
  //         (morton_id == b.morton_id && id < b.id);
  // }

  /* return Z value of this point */
  unsigned long long interleave_bits() const {
    // Pun the x and y coordinates as integers: Just re-interpret the bits.
    //
    auto ix = static_cast<unsigned int>(x);
    auto iy = static_cast<unsigned int>(y);
    // cout << ix << ", " << iy << endl;
    // cout << bitset<32>(ix) << endl;
    // cout << bitset<32>(iy) << endl;

    auto ret = 0ull;
    for (auto i = 0; i < 32; i++) {
      ret |= ((ix & (1ull << i)) << (i + 1)) | ((iy & (1ull << i)) << i);
    }
    // cout << bitset<64>(ret) << endl;
    return ret;
  }

  unsigned long long overlap_bits() const {
    auto ix = static_cast<unsigned long long>(x);
    auto iy = static_cast<unsigned long long>(y);
    unsigned long long p[] = {ix, iy};
    return hilbert_c2i(2, 32, p);
  }

  long long mortonIndex() const {
    // Pun the x and y coordinates as integers: Just re-interpret the bits.
    //
    auto ix = static_cast<unsigned int>(x);
    auto iy = static_cast<unsigned int>(y);
    // cout << ix << " " << iy << endl;

    // Since we're assuming 2s complement arithmetic (99.99% of hardware today),
    // we'll need to convert these raw integer-punned floats into
    // their corresponding integer "indices".

    // Smear their sign bits into these for twiddling below.
    //
    auto const ixs = static_cast<int>(ix) >> 31;
    auto const iys = static_cast<int>(iy) >> 31;

    // This is a combination of a fast absolute value and a bias.
    //
    // We need to adjust the values so -FLT_MAX is close to 0.
    //
    ix = (((ix & 0x7FFFFFFFL) ^ ixs) - ixs) + 0x7FFFFFFFL;
    iy = (((iy & 0x7FFFFFFFL) ^ iys) - iys) + 0x7FFFFFFFL;

    // Now we have -FLT_MAX close to 0, and FLT_MAX close to UINT_MAX,
    // with everything else in-between.
    //
    // To make this easy, we'll work with x and y as 64-bit integers.
    //
    long long xx = ix;
    long long yy = iy;

    // Dilate and combine as usual...

    xx = (xx | (xx << 16)) & 0x0000ffff0000ffffLL;
    yy = (yy | (yy << 16)) & 0x0000ffff0000ffffLL;

    xx = (xx | (xx << 8)) & 0x00ff00ff00ff00ffLL;
    yy = (yy | (yy << 8)) & 0x00ff00ff00ff00ffLL;

    xx = (xx | (xx << 4)) & 0x0f0f0f0f0f0f0f0fLL;
    yy = (yy | (yy << 4)) & 0x0f0f0f0f0f0f0f0fLL;

    xx = (xx | (xx << 2)) & 0x3333333333333333LL;
    yy = (yy | (yy << 2)) & 0x3333333333333333LL;

    xx = (xx | (xx << 1)) & 0x5555555555555555LL;
    yy = (yy | (yy << 1)) & 0x5555555555555555LL;

    return xx | (yy << 1);
  }
};

typedef pair<Point, Point> Bounding_Box;

template <class T>
auto read_pts(T& P, ifstream& fin, bool real_data = false) {
  if (!real_data) {
    size_t n, d;
    fin >> n >> d;
    P.resize(n);
    size_t id;
    FT x, y;
    FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX),
        y_max(FT_INF_MIN);
    for (size_t i = 0; i < n; i++) {
      // fin >> id >> x >> y;
      fin >> x >> y;
      x_max = max(x_max, x);
      x_min = min(x_min, x);
      y_max = max(y_max, y);
      y_min = min(y_min, y);
      id = i;
      auto cur_p = Point(id, x, y);
      P[i] = cur_p;
    }
    return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
  } else {
    size_t id;
    FT x, y;
    FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX),
        y_max(FT_INF_MIN);
    P.clear();
    while (fin >> id >> x >> y) {
      x *= 1000000;
      y *= 1000000;
      if (x < 0 || y < 0) {
        // cout << "negative points found" << endl;
        // cout << x << " " << y << endl;
        // exit(-1);
      }
      x_max = max(x_max, x);
      x_min = min(x_min, x);
      y_max = max(y_max, y);
      y_min = min(y_min, y);
      auto cur_p = Point(id, x, y);
      P.emplace_back(cur_p);
    }
    return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
  }
}

template <class T>
void print_binary(T x) {
  cout << bitset<sizeof(x) * 8>(x) << endl;
}

template <class MBR>
void print_mbr(MBR& mbr) {
  cout << fixed << setprecision(6) << "[(" << mbr.first.x << ", " << mbr.first.y
       << "), (" << mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
}

template <class MBR>
bool is_same_mbr(MBR& lhs, MBR& rhs) {
  // cout << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" <<
  // mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
  return lhs.first.x == rhs.first.x && lhs.first.y == rhs.first.y &&
         lhs.second.x == rhs.second.x && lhs.second.y == rhs.second.y;
}

template <class T>
size_t split_by_bit(T& P, size_t l, size_t r, size_t b) {
  // size_t split_value = (1ull << (b - 1));
  unsigned int splitter = (1u << ((b - 1) / 2));
  auto less = [&](auto pt) {
    return b % 2 ? (static_cast<unsigned int>(pt.y) & splitter) == 0
                 : (static_cast<unsigned int>(pt.x) & splitter) == 0;
  };
  size_t start = l, end = r, mid = start;
  while (end - start > 16) {
    mid = (start + end) / 2;
    if (!less(P[mid]))
      end = mid;
    else
      start = mid + 1;
  }
  mid = end;
  for (auto i = start; i < end; i++) {
    if (!less(P[i])) {
      mid = i;
      return i;
    }
  }
  return mid;
}

// merge two bounding boxes, please make sure *both bounding boxes are valid*.
template <class MBR>
MBR merge_mbr(MBR a, MBR b) {
  return {Point(min(a.first.x, b.first.x), min(a.first.y, b.first.y)),
          Point(max(a.second.x, b.second.x), max(a.second.y, b.second.y))};
}

// check whether a given point inside a mbr (boundary included).
template <class MBR>
bool point_in_mbr(Point p, MBR mbr) {
  return mbr.first.x <= p.x && p.x <= mbr.second.x && mbr.first.y <= p.y &&
         p.y <= mbr.second.y;
}

// check whether a given mbr inside a mbr, e,g., check whether small_mbr in
// large_mbr (boundary included).
template <class MBR>
bool mbr_in_mbr(MBR& small_mbr, MBR& large_mbr) {
  return small_mbr.first.x >= large_mbr.first.x &&
         small_mbr.second.x <= large_mbr.second.x &&
         small_mbr.first.y >= large_mbr.first.y &&
         small_mbr.second.y <= large_mbr.second.y;
}

// check whether two given mbrs exclude each other (do not check boundaries).
template <class MBR>
bool mbr_exclude_mbr(MBR& small_mbr, MBR& large_mbr) {
  return (small_mbr.first.x > large_mbr.second.x ||
          small_mbr.second.x < large_mbr.first.x) ||
         (small_mbr.first.y > large_mbr.second.y ||
          small_mbr.second.y < large_mbr.first.y);
}

// relations between two mbrs: -1: excluded; 0: intersected; 1: the smaller one
// is contained by the larger one;
template <class MBR>
int mbr_mbr_relation(MBR& small_mbr, MBR& large_mbr) {
  auto minc_x = max(small_mbr.first.x, large_mbr.first.x);
  auto minc_y = max(small_mbr.first.y, large_mbr.first.y);
  auto maxc_x = min(small_mbr.second.x, large_mbr.second.x);
  auto maxc_y = min(small_mbr.second.y, large_mbr.second.y);
  if (minc_x <= maxc_x && minc_y <= maxc_y) {
    if (minc_x == small_mbr.first.x && maxc_x == small_mbr.second.x &&
        minc_y == small_mbr.first.y && maxc_y == small_mbr.second.y)
      return 1;
    return 0;
  }
  return -1;
}

template <class Records>
auto get_mbr(Records& P) {
  FT x_min = FT_INF_MAX, x_max = FT_INF_MIN, y_min = FT_INF_MAX,
     y_max = FT_INF_MIN;
  for (auto& p : P) {
    x_min = min(x_min, p.x);
    x_max = max(x_max, p.x);
    y_min = min(y_min, p.y);
    y_max = max(y_max, p.y);
  }
  return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
}

// return the sqr distance between a point and a mbr
template <class MBR>
auto point_mbr_sqrdis(Point p, MBR& mbr) {
  FT dx =
      max(max(mbr.first.x - p.x, (FT)0.0), max(p.x - mbr.second.x, (FT)0.0));
  FT dy =
      max(max(mbr.first.y - p.y, (FT)0.0), max(p.y - mbr.second.y, (FT)0.0));
  return dx * dx + dy * dy;
}

// return the sqr distance between two points
auto point_point_sqrdis(Point& lhs, Point& rhs) {
  return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y);
}

// check whether two given MBRs are intersected.
template <class MBR>
bool mbr_intersects_mbr(MBR a, MBR b) {
  if (a.first.x > b.second.x || a.second.x < b.first.x) return false;
  if (a.first.y > b.second.y || a.second.y < b.first.y) return false;
  return true;
}

template <class T>
auto generate_range_query(T& P, size_t n) {
  auto id1 = rand() % n;
  auto id2 = rand() % n;
  FT xmin = min(P[id1].x, P[id2].x), xmax = max(P[id1].x, P[id2].x),
     ymin = min(P[id1].y, P[id2].y), ymax = max(P[id1].y, P[id2].y);
  return make_pair(Point(xmin, ymin), Point(xmax, ymax));
}

template <class In>
auto read_range_query(In qry_in) {
  ifstream fin(qry_in);
  size_t n, d;
  fin >> n >> d;
  parlay::sequence<Bounding_Box> ret(n);
  for (size_t i = 0; i < n; i++) {
    fin >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >>
        ret[i].second.y;
  }
  return ret;
}

template <class In>
auto read_range_query(In qry_in, size_t q_type, size_t& maxSize) {
  ifstream fin(qry_in);
  if (q_type == 8) {  // range report, need maxSize
    fin >> maxSize;
  }
  size_t n, d;
  fin >> n >> d;
  parlay::sequence<Bounding_Box> ret(n);
  for (size_t i = 0; i < n; i++) {
    fin >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >>
        ret[i].second.y;
  }
  return ret;
}

// data type/helper functions for nearest neighbor search
typedef pair<Point, FT> nn_pair;

struct nn_pair_cmp {
  bool operator()(nn_pair& lhs, nn_pair& rhs) {
    return lhs.second < rhs.second ||
           (lhs.second == rhs.second && lhs.first.id > rhs.first.id);
  }
};

template <class Pset>
FT knn_bf(size_t& k, Point q, Pset& P) {
  vector<FT> q_sqrdis = {};
  for (size_t i = 0; i < P.size(); i++) {
    auto cur_sqrdis = point_point_sqrdis(q, P[i]);
    q_sqrdis.emplace_back(cur_sqrdis);
  }
  sort(q_sqrdis.begin(), q_sqrdis.end());
  return q_sqrdis[k - 1];
}

template <class Pset>
void morton_sort(Pset& P) {
  sort(P.begin(), P.end(), [&](auto lhs, auto rhs) {
    auto msd = 0;
    if (geobase::less_msb(
            static_cast<unsigned int>(lhs.x) ^ static_cast<unsigned int>(rhs.x),
            static_cast<unsigned int>(lhs.y) ^
                static_cast<unsigned int>(rhs.y)))
      msd = 1;
    return !msd ? lhs.x < rhs.x : lhs.y < rhs.y;
  });
}

template <class P>
bool morton_less(P const& lhs, P const& rhs) {
  auto msd = 0;
  if (geobase::less_msb(
          static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x),
          static_cast<unsigned int>(lhs->y) ^
              static_cast<unsigned int>(rhs->y)))
    msd = 1;
  return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
}

template <typename Pset>
auto morton_merge(Pset& lhs, Pset& rhs) {
  size_t i = 0, j = 0, k = 0, n = lhs.size(), m = rhs.size();
  Pset ret;
  ret.resize(n + m);
  while (i < n && j < m) {
    if (morton_less(lhs[i], rhs[j])) {
      ret[k++] = lhs[i++];
    } else {
      ret[k++] = rhs[j++];
    }
  }
  while (i < n) ret[k++] = lhs[i++];
  while (j < m) ret[k++] = rhs[j++];
  return ret;
}

// generate n random points within the largest bounding box
template <typename Pset>
auto random_sample(size_t n, Pset P) {
  unsigned seed = 233666;
  default_random_engine e(seed);
  vector<size_t> sampled_ids(P.size());
  parlay::parallel_for(0, P.size(), [&](size_t i) { sampled_ids[i] = i; });
  shuffle(sampled_ids.begin(), sampled_ids.end(), e);
  auto P2 = P.substr(0, n);
  for (size_t i = 0; i < n; i++) {
    P2[i] = P[sampled_ids[i]];
  }
  return P2;
}

//  return a set of sorted points
template <typename PT>
auto get_sorted_points(PT& P, bool use_hilbert = false) {
  auto n = P.size();
  parlay::parallel_for(0, n, [&](int i) {
    P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
  });

  return parlay::sort(
      P, [&](auto lhs, auto rhs) { return lhs.morton_id < rhs.morton_id; });
  // if (use_hilbert){
  // }
  // else{
  //     return parlay::sort(P,  [&](auto lhs, auto rhs){
  // 	    auto msd = 0;
  // 	    if (geobase::less_msb(static_cast<unsigned int>(lhs.x) ^
  // static_cast<unsigned int>(rhs.x), static_cast<unsigned int>(lhs.y) ^
  // static_cast<unsigned int>(rhs.y))) 		    msd = 1;
  // return !msd ? lhs.x < rhs.x : lhs.y < rhs.y;
  //     });
  // }
}

template <typename PT>
auto get_sorted_points_hilbert(PT& P, bool use_hilbert = true) {
  auto n = P.size();
  parlay::parallel_for(0, n, [&](int i) {
    P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
  });

  return parlay::sort(
      P, [&](auto lhs, auto rhs) { return lhs.morton_id < rhs.morton_id; });
  // auto P_set = parlay::sort(P, [&](auto lhs, auto rhs){
  // 	unsigned long long p_lhs[] = {static_cast<unsigned long long>(lhs.x),
  // static_cast<unsigned long long>(lhs.y)}; 	unsigned long long p_rhs[] =
  // {static_cast<unsigned long long>(rhs.x), static_cast<unsigned long
  // long>(rhs.y)}; 	return hilbert_cmp(2, sizeof(unsigned long long), 8 *
  // sizeof(unsigned long long), p_lhs, p_rhs) == -1;
  // });
  // return P_set;
}

//	return a set of sorted address, do not modify the input;
template <typename PT>
auto get_sorted_address(PT& P) {
  parlay::sequence<geobase::Point*> P_set(P.size());
  parlay::parallel_for(0, P.size(), [&](int i) { P_set[i] = &P[i]; });
  P_set = parlay::sort(P_set, [&](auto lhs, auto rhs) {
    auto msd = 0;
    if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^
                              static_cast<unsigned int>(rhs->x),
                          static_cast<unsigned int>(lhs->y) ^
                              static_cast<unsigned int>(rhs->y)))
      msd = 1;
    return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
  });
  return P_set;
}

template <typename PT>
auto get_unsorted_address(PT& P) {
  parlay::sequence<Point*> P_set(P.size());
  parlay::parallel_for(0, P.size(), [&](int i) { P_set[i] = &P[i]; });
  return P_set;
}

template <typename PT>
auto record_check(PT& P) {
  map<int, int> mmp = {};
  auto cnt = 0;
  auto num = 0;
  for (size_t i = 0; i < P.size(); i++) {
    mmp[P[i]->id]++;
  }
  for (auto& key_val : mmp) {
    if (key_val.second > 1) {
      num++;
      cnt += key_val.second;
      // cout << key_val.first << ", " << key_val.second << endl;
    }
  }
  return make_tuple(num, cnt);
}

// return a set contains all ids in a sequence of point pointer
template <typename PT>
auto get_point_id(PT& P) {
  set<int> ret = {};
  for (auto pt : P) {
    ret.insert(pt->id);
  }
  return ret;
}

template <typename PT>
auto count_duplicate_zvalues(PT& P) {
  unordered_map<unsigned long long, int> zvalue_map = {};
  for (auto& pt : P) {
    zvalue_map[pt.morton_id]++;
  }
  auto cnt = 0, total_cnt = 0;
  for (auto const& pair : zvalue_map) {
    if (pair.second > 1) {
      // cout << "Number " << pair.first << " appears " << pair.second << "
      // times." << endl;
      cnt++;
      total_cnt += pair.second;
    }
  }
  cout << "Total # of duplicate keys: " << cnt << ", " << total_cnt << endl;
  return cnt;
}

//  merge to sequence by id, return two sequences. The first one has no
//  conflict, the second one has conflict
template <typename PT>
auto merge_by_id_with_conflict(PT& a, PT& b) {
  size_t i = 0, j = 0;
  auto id_cmp = [&](auto lhs, auto rhs) { return lhs.id < rhs.id; };
  auto sorted_a = parlay::sort(a, id_cmp);
  auto sorted_b = parlay::sort(b, id_cmp);
  PT ret_noconflict = {}, ret_conflict = {};
  while (i < sorted_a.size() && j < sorted_b.size()) {
    if (sorted_a[i].id == sorted_b[j].id) {  //  two points have the same id
      if (!(sorted_a[i] ==
            sorted_b[j])) {  //  coordinate does not match, conflict
        ret_conflict.emplace_back(sorted_a[i++]);
        ret_conflict.emplace_back(sorted_b[j++]);
      } else {  //  no conflict
        ret_noconflict.emplace_back(sorted_a[i++]);
        j++;
      }
    } else {  //  no conflict, store the one with smaller id
      if (sorted_a[i].id < sorted_b[j].id) {
        ret_noconflict.emplace_back(sorted_a[i++]);
      } else {
        ret_noconflict.emplace_back(sorted_b[j++]);
      }
    }
  }
  while (i < sorted_a.size()) ret_noconflict.emplace_back(sorted_a[i++]);
  while (j < sorted_b.size()) ret_noconflict.emplace_back(sorted_b[j++]);
  return make_tuple(ret_noconflict, ret_conflict);
}
}  // namespace geobase
