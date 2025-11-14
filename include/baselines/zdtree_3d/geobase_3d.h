#pragma once

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

#include "parlay/internal/merge.h"
#include "psi/dependence/comparator.h"
#include "psi/dependence/space_filling_curve/hilbert.h"
#include "psi/dependence/splitter.h"

namespace ZD3D {
namespace geobase {
using namespace std;
using namespace psi::hilbert;
using FT = double;
// using FT = float;
constexpr FT FT_INF_MIN = numeric_limits<FT>::min();
constexpr FT FT_INF_MAX = numeric_limits<FT>::max();
constexpr FT FT_EPS = numeric_limits<FT>::epsilon();

inline int dcmp(const FT& x) {
  if (fabs(x) < FT_EPS) return 0;
  return x < 0 ? -1 : 1;
}

bool less_msb(unsigned int x, unsigned int y) { return x < y && x < (x ^ y); }

struct Point {
  size_t id = 0;
  FT x, y, z;  // add z
  unsigned long long morton_id;
  Point() {}
  Point(FT _x, FT _y) : x(_x), y(_y) {}
  Point(size_t _id, FT _x, FT _y) : id(_id), x(_x), y(_y) {}

  // 3D
  Point(FT _x, FT _y, FT _z) : x(_x), y(_y), z(_z) {}
  Point(size_t _id, FT _x, FT _y, FT _z) : id(_id), x(_x), y(_y), z(_z) {}

  // new things
  using Coord = FT;
  // using Coords = std::array<Coord, 2>;
  using Coords = std::array<Coord, 3>;
  using Num = psi::Num_Comparator<Coord>;
  using DisType =
      std::conditional_t<std::is_integral_v<Coord>, int_fast64_t, double>;
  using DimsType = uint_fast8_t;
  using BP = Point;
  Point(Coords const& _pnt) : x(_pnt[0]), y(_pnt[1]), z(_pnt[2]) {}
  Point(FT x) : x(x), y(x), z(x) {}
  Coord& operator[](DimsType i) {
    if (i == 0) return x;
    return i == 1 ? y : z;
  }
  Coord const& operator[](DimsType i) const {
    if (i == 0) return x;
    return i == 1 ? y : z;
  }
  Point const MinCoords(Point const& b) const {
    Point p;
    p[0] = Num::Min(x, b[0]);
    p[1] = Num::Min(y, b[1]);
    p[2] = Num::Min(z, b[2]);
    return p;
  }

  Point const MaxCoords(Point const& b) const {
    Point p;
    p[0] = Num::Max(x, b[0]);
    p[1] = Num::Max(y, b[1]);
    p[2] = Num::Max(z, b[2]);
    return p;
  }
  // static consteval auto GetDim() { return 2; }
  static consteval auto GetDim() { return 3; }
  bool SameDimension(Point const& b) const {
    return x == b.x && y == b.y && z == b.z;
  }
  Coords GetCoords() const { return {x, y, z}; }

  // original things
  bool operator==(Point const& p) const {
    return !(id - p.id) && !dcmp(x - p.x) && !dcmp(y - p.y) && !dcmp(z - p.z);
  }

  friend std::ostream& operator<<(std::ostream& os, Point const& p) {
    // os << fixed << setprecision(6) << p.id << ": (" << p.x << ", " << p.y <<
    // ", " << p.z << ")" << bitset<64>(p.morton_id); os << fixed <<
    // setprecision(6) << p.id << ": (" << p.x << ", " << p.y << ", " << p.z <<
    // ")";
    os << fixed << setprecision(6) << "Zvalue: " << bitset<64>(p.morton_id)
       << ", ID: " << p.id << ", XYZ: (" << p.x << ", " << p.y << ", " << p.z
       << ")";
    return os;
  }

  /* return Z value of this point */
  unsigned long long interleave_bits() const {
    /* 3D Morton Value Form: xyzxyz...xyz- the lowest bit is not used. */
    auto ix = static_cast<unsigned int>(x);
    auto iy = static_cast<unsigned int>(y);
    auto iz = static_cast<unsigned int>(z);
    unsigned long long ret = 0ull;
    // std::cout << std::bitset<32>(ix) << std::endl;
    // std::cout << std::bitset<32>(iy) << std::endl;
    // std::cout << std::bitset<32>(iz) << std::endl;
    for (int i = 0; i < 21; i++) {
      // std::cout << "i = " << i << std::endl;
      // ret |= (((ix & (1ull << i)) << (3 * i + 2)) | ((iy & (1ull << i)) << (3
      // * i + 1)) | ((iz & (1ull << i)) << (3 * i)));
      ret |= ((iz >> i) & 1ull) << (3 * i);
      ret |= ((iy >> i) & 1ull) << (3 * i + 1);
      ret |= ((ix >> i) & 1ull) << (3 * i + 2);
      // std::cout << std::bitset<64>(ret) << std::endl;
    }

    // ret <<= 1;
    return ret;

    // /* Original Implementation for 2D */
    // // Pun the x and y coordinates as integers: Just re-interpret the bits.
    // auto ix = static_cast<unsigned int>(x);
    // auto iy = static_cast<unsigned int>(y);

    // auto ret = 0ull;
    // for (auto i = 0; i < 32; i++) {
    //   ret |= ((ix & (1ull << i)) << (i + 1)) | ((iy & (1ull << i)) << i);
    // }
    // return ret;
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

struct diff_type {
  size_t add_cnt, remove_cnt;
  parlay::sequence<Point> add, remove;
  diff_type() {
    add_cnt = 0;
    remove_cnt = 0;
  }
  diff_type(size_t add_sz, size_t remove_sz) {
    add_cnt = 0;
    remove_cnt = 0;
    add = parlay::sequence<Point>::uninitialized(add_sz);
    remove = parlay::sequence<Point>::uninitialized(remove_sz);
  }
  void add_point(Point& p, bool reverse = false) {
    if (!reverse) {
      parlay::assign_uninitialized(add[add_cnt++], p);
    } else {
      parlay::assign_uninitialized(remove[remove_cnt++], p);
    }
  }
  void remove_point(Point& p, bool reverse = false) {
    if (!reverse) {
      parlay::assign_uninitialized(remove[remove_cnt++], p);
    } else {
      parlay::assign_uninitialized(add[add_cnt++], p);
    }
  }
  void compact() {
    add.resize(add_cnt);
    remove.resize(remove_cnt);
  }
  void reset() {
    add_cnt = 0;
    remove_cnt = 0;
  }
  void reset(size_t add_sz, size_t remove_sz) {
    add_cnt = 0;
    remove_cnt = 0;
    add = parlay::sequence<Point>::uninitialized(add_sz);
    remove = parlay::sequence<Point>::uninitialized(remove_sz);
  }
};

// template <class T>
// auto read_pts(T& P, ifstream& fin, bool real_data = false) {
//   if (!real_data) {
//     size_t n, d;
//     fin >> n >> d;
//     P.resize(n);
//     size_t id;
//     FT x, y;
//     FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX),
//         y_max(FT_INF_MIN);
//     for (size_t i = 0; i < n; i++) {
//       // fin >> id >> x >> y;
//       fin >> x >> y;
//       x_max = max(x_max, x);
//       x_min = min(x_min, x);
//       y_max = max(y_max, y);
//       y_min = min(y_min, y);
//       id = i;
//       auto cur_p = Point(id, x, y);
//       P[i] = cur_p;
//     }
//     return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
//   } else {
//     size_t id;
//     FT x, y;
//     FT x_min(FT_INF_MAX), x_max(FT_INF_MIN), y_min(FT_INF_MAX),
//         y_max(FT_INF_MIN);
//     P.clear();
//     while (fin >> id >> x >> y) {
//       x *= 1000000;
//       y *= 1000000;
//       if (x < 0 || y < 0) {
//         continue;  // ignore outliers
//       }
//       x_max = max(x_max, x);
//       x_min = min(x_min, x);
//       y_max = max(y_max, y);
//       y_min = min(y_min, y);
//       auto cur_p = Point(id, x, y);
//       P.emplace_back(cur_p);
//     }
//     return Bounding_Box({Point(x_min, y_min), Point(x_max, y_max)});
//   }
// }

template <typename PT>
auto filter_diff_results(PT& add, PT& remove) {
  auto id_cmp = [&](auto lhs, auto rhs) { return lhs.id < rhs.id; };
  auto sorted_add = parlay::sort(add, id_cmp);
  auto sorted_remove = parlay::sort(remove, id_cmp);

  PT insert_points = {}, delete_points = {}, update_points = {};
  size_t i = 0, j = 0;
  while (i < sorted_add.size() && j < sorted_remove.size()) {
    if (sorted_add[i].id <
        sorted_remove[j].id) {  //	we should insert sorted_add[i]
      insert_points.emplace_back(sorted_add[i++]);
    } else if (sorted_add[i].id ==
               sorted_remove[j].id) {  //	exist in both add and remove,
                                       // indicate it is an update point
      update_points.emplace_back(sorted_add[i++]);
      j++;
    } else {  //	we should remove sorted_remove[j]
      delete_points.emplace_back(sorted_remove[j++]);
    }
  }
  while (i < sorted_add.size()) insert_points.emplace_back(sorted_add[i++]);
  while (j < sorted_remove.size())
    delete_points.emplace_back(sorted_remove[j++]);
  return make_tuple(insert_points, delete_points, update_points);
}

template <class T>
void print_binary(T x) {
  cout << bitset<sizeof(x) * 8>(x) << endl;
}

template <class MBR>
void print_mbr(MBR& mbr) {
  cout << fixed << setprecision(6) << "[(" << mbr.first.x << ", " << mbr.first.y
       << ", " << mbr.first.z << "), (" << mbr.second.x << ", " << mbr.second.y
       << ", " << mbr.second.z << ")" << "]" << endl;
}

template <class MBR>
bool is_same_mbr(MBR& lhs, MBR& rhs) {
  // cout << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" <<
  // mbr.second.x << ", " << mbr.second.y << ")" << "]" << endl;
  return lhs.first.x == rhs.first.x && lhs.first.y == rhs.first.y &&
         lhs.first.z == rhs.first.z && lhs.second.x == rhs.second.x &&
         lhs.second.y == rhs.second.y && lhs.second.z == rhs.second.z;
}

template <class T>
size_t split_by_bit(T& P, size_t l, size_t r, size_t b) {
  unsigned int splitter = (b - 1) / 3;
  // auto less = [&](auto pt) {
  //   return b % 2 ? (static_cast<unsigned int>(pt.y) & splitter) == 0
  //                : (static_cast<unsigned int>(pt.x) & splitter) == 0;
  // };
  auto less = [&](auto pt) {
    if (b % 3 == 0) {  // x axis
      auto x = static_cast<unsigned int>(pt.x);
      return ((x >> splitter) & 1) == 0;
    } else if (b % 3 == 2) {  // y axis
      auto y = static_cast<unsigned int>(pt.y);
      return ((y >> splitter) & 1) == 0;
    } else {
      auto z = static_cast<unsigned int>(pt.z);
      return ((z >> splitter) & 1) == 0;
    }
  };
  size_t start = l, end = r, mid = start;
  while (end - start > 16) {
    mid = (start + end) / 2;
    if (less(P[mid]))
      start = mid + 1;
    else
      end = mid;
  }
  mid = end;
  // cout << "start = " << start << ", end = " << end << endl;
  for (auto i = start; i < end; i++) {
    // cout << "i = " << i << endl;
    // if (b == 4){
    //   auto z = static_cast<unsigned int>(P[i].z);
    //   cout << "z value = " << z << endl;
    //   cout << ((z >> 1) & 1) << endl;
    // }
    if (!less(P[i])) {
      mid = i;
      return i;
    }
  }
  return mid;
}

// merge two bounding boxes, please make sure *both bounding boxes are valid*.
template <class MBR>
MBR merge_mbr(MBR& a, MBR& b) {
  return {Point(min(a.first.x, b.first.x), min(a.first.y, b.first.y),
                min(a.first.z, b.first.z)),
          Point(max(a.second.x, b.second.x), max(a.second.y, b.second.y),
                max(a.second.z, b.second.z))};
}

// check whether a given point inside a mbr (boundary included).
template <class MBR>
bool point_in_mbr(Point& p, MBR& mbr) {
  return mbr.first.x <= p.x && p.x <= mbr.second.x && mbr.first.y <= p.y &&
         p.y <= mbr.second.y && mbr.first.z <= p.z && p.z <= mbr.second.z;
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
int mbr_mbr_relation(MBR& small_mbr, MBR& large_mbr, size_t d = 3) {
  auto minc_x = max(small_mbr.first.x, large_mbr.first.x);
  auto maxc_x = min(small_mbr.second.x, large_mbr.second.x);
  auto minc_y = max(small_mbr.first.y, large_mbr.first.y);
  auto maxc_y = min(small_mbr.second.y, large_mbr.second.y);
  if (d == 2) {
    if (minc_x <= maxc_x && minc_y <= maxc_y) {
      if (minc_x == small_mbr.first.x && maxc_x == small_mbr.second.x &&
          minc_y == small_mbr.first.y && maxc_y == small_mbr.second.y)
        return 1;
      return 0;
    }
    return -1;
  }

  /* 3D case */
  auto minc_z = max(small_mbr.first.z, large_mbr.first.z);
  auto maxc_z = min(small_mbr.second.z, large_mbr.second.z);
  if (minc_x <= maxc_x && minc_y <= maxc_y && minc_z <= maxc_z) {
    if (minc_x == small_mbr.first.x && maxc_x == small_mbr.second.x &&
        minc_y == small_mbr.first.y && maxc_y == small_mbr.second.y &&
        minc_z == small_mbr.first.z && maxc_z == small_mbr.second.z)
      return 1;
    return 0;
  }
  return -1;
}

template <class Records>
auto get_mbr(Records& P) {
  FT x_min = FT_INF_MAX, x_max = FT_INF_MIN, y_min = FT_INF_MAX,
     y_max = FT_INF_MIN, z_min = FT_INF_MAX, z_max = FT_INF_MIN;
  for (auto& p : P) {
    x_min = min(x_min, p.x);
    x_max = max(x_max, p.x);
    y_min = min(y_min, p.y);
    y_max = max(y_max, p.y);
    z_min = min(z_min, p.z);
    z_max = max(z_max, p.z);
  }
  return Bounding_Box({Point(x_min, y_min, z_min), Point(x_max, y_max, z_max)});
}

auto mbr_mbr_within_dis(Bounding_Box& mbr1, Bounding_Box& mbr2, FT& point_dis) {
  return !(dcmp(mbr1.second.x + point_dis - mbr2.first.x) < 0 ||
           dcmp(mbr1.second.y + point_dis - mbr2.first.y) < 0 ||
           dcmp(mbr2.second.x + point_dis - mbr1.first.x) < 0 ||
           dcmp(mbr2.second.y + point_dis - mbr1.first.y) < 0);
}

// return the sqr distance between a point and a mbr
template <class MBR>
auto point_mbr_sqrdis(Point& p, MBR& mbr) {
  FT dx =
      max(max(mbr.first.x - p.x, (FT)0.0), max(p.x - mbr.second.x, (FT)0.0));
  FT dy =
      max(max(mbr.first.y - p.y, (FT)0.0), max(p.y - mbr.second.y, (FT)0.0));
  FT dz =
      max(max(mbr.first.z - p.z, (FT)0.0), max(p.z - mbr.second.z, (FT)0.0));

  return dx * dx + dy * dy + dz * dz;
}

// return the sqr distance between two points
auto point_point_sqrdis(Point& lhs, Point& rhs) {
  return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y) +
         (lhs.z - rhs.z) * (lhs.z - rhs.z);
}

/* Add Ziyang's box check logic */
template <class MBR>
bool box_same_box(MBR& a, MBR& b) {
  return !dcmp(a.first.x - b.first.x) && !dcmp(a.first.y - b.first.y) &&
         !dcmp(a.first.z - b.first.z) && !dcmp(a.second.x - b.second.x) &&
         !dcmp(a.second.y - b.second.y) && !dcmp(a.second.z - b.second.z);
}

template <class MBR>
bool box_within_box(MBR& a, MBR& b) {
  return dcmp(a.first.x - b.first.x) >= 0 && dcmp(a.first.y - b.first.y) >= 0 &&
         dcmp(a.first.z - b.first.z) >= 0 &&
         dcmp(a.second.x - b.second.x) <= 0 &&
         dcmp(a.second.y - b.second.y) <= 0 &&
         dcmp(a.second.z - b.second.z) <= 0;
}

template <class MBR>
bool box_intersect_box(MBR& a, MBR& b) {
  return !(
      dcmp(a.second.x - b.first.x) < 0 || dcmp(a.first.x - b.second.x) > 0 ||
      dcmp(a.second.y - b.first.y) < 0 || dcmp(a.first.y - b.second.y) > 0 ||
      dcmp(a.second.z - b.first.z) < 0 || dcmp(a.first.z - b.second.z) > 0);
}

template <class MBR>
auto zy_mbr_relation(MBR& a, MBR& b) {
  if (!box_intersect_box(a, b)) return -1;
  if (box_within_box(a, b)) return 1;
  return 0;
}

// check whether two given MBRs are intersected.
template <class MBR>
bool mbr_intersects_mbr(MBR& a, MBR& b) {
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

// template <class In>
// auto read_range_query(In qry_in) {
//   ifstream fin(qry_in);
//   size_t n, d;
//   fin >> n >> d;
//   parlay::sequence<Bounding_Box> ret(n);
//   for (size_t i = 0; i < n; i++) {
//     fin >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >>
//         ret[i].second.y;
//   }
//   return ret;
// }

// template <class In>
// auto read_range_query(In qry_in, size_t q_type, size_t& maxSize,
//                       bool is_real = false) {
//   ifstream fin(qry_in);
//   if (q_type == 8) {  // range report, need maxSize
//     fin >> maxSize;
//   }
//   size_t n, d;
//   fin >> n >> d;
//   parlay::sequence<Bounding_Box> ret(n);
//   parlay::sequence<size_t> cnt(n);
//   for (size_t i = 0; i < n; i++) {
//     fin >> cnt[i] >> ret[i].first.x >> ret[i].first.y >> ret[i].second.x >>
//         ret[i].second.y;
//   }
//   return make_tuple(cnt, ret);
// }

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

template <typename Pset>
auto shuffle_point(Pset& P, size_t substr_size = 0) {
  auto n = P.size();
  parlay::sequence<Point> rand_p(n);
  // random_device rd;
  auto rand_seed = 233666;
  mt19937 g(rand_seed);
  parlay::parallel_for(0, n, [&](size_t i) { rand_p[i] = P[i]; });
  shuffle(rand_p.begin(), rand_p.end(), g);
  if (substr_size > 0) {
    rand_p = rand_p.substr(0, substr_size);
  }
  return rand_p;
}

/* Insert ratio is from 0 to 10. Insert point will add bias to id to make them
 * unique. */
template <typename Pset>
auto split_insert_delete(Pset& P, size_t& insert_ratio, size_t id_bias) {
  size_t insert_num = P.size() / 10 * insert_ratio;
  size_t delete_num = P.size() - insert_num;
  auto P_insert = P.substr(0, insert_num);
  parlay::parallel_for(0, insert_num,
                       [&](size_t i) { P_insert[i].id += id_bias; });
  auto P_delete = P.substr(P.size() - delete_num, delete_num);
  return make_tuple(P_insert, P_delete);
}

template <typename Pset>
auto print_Pset_info(Pset& P, string name = "debug", size_t end_idx = 0) {
  cout << name << ": " << endl;
  if (end_idx == 0) end_idx = P.size();
  for (size_t i = 0; i < end_idx; i++) {
    cout << P[i] << endl;
  }
}

template <typename Pset>
auto collect_newver_point(Pset& P, Pset& P_insert, Pset& P_delete) {
  unordered_map<size_t, bool> removed = {};
  for (auto& pt : P_delete) {
    removed[pt.id] = true;
  }
  auto ret = parlay::sequence<Point>::uninitialized(P.size() + P_insert.size());
  auto cnt = 0;
  for (auto& pt : P) {
    if (!removed[pt.id]) parlay::assign_uninitialized(ret[cnt++], pt);
  }
  for (auto& pt : P_insert) {
    // if (!removed[pt.id])
    parlay::assign_uninitialized(ret[cnt++], pt);
  }
  ret.resize(cnt);
  return ret;
}

auto compute_cur_box(Bounding_Box& cur_mbr, FT& x_prefix, FT& y_prefix,
                     size_t b, bool x_splitter) {
  size_t shift_b = (b + 1) / 2;
  FT split_value = 1.0 * (1u << (shift_b - 1));

  auto L_box = cur_mbr;
  auto R_box = cur_mbr;
  auto rx_prefix = x_prefix;
  auto ry_prefix = y_prefix;
  if (x_splitter) {
    L_box.second.x = min(x_prefix + split_value - FT_EPS, L_box.second.x);
    R_box.first.x = max(x_prefix + split_value, R_box.first.x);
    rx_prefix += split_value;
  } else {
    L_box.second.y = min(y_prefix + split_value - FT_EPS, L_box.second.y);
    R_box.first.y = max(y_prefix + split_value, R_box.first.y);
    ry_prefix += split_value;
  }
  return make_tuple(L_box, R_box, rx_prefix, ry_prefix);
}

auto compute_cur_box3(Bounding_Box& cur_mbr, FT& x_prefix, FT& y_prefix,
                      FT& z_prefix, size_t b, size_t split_axis) {
  size_t shift_b = (b + 1) / 3;
  FT split_value = 1.0 * (1u << (shift_b - 1));
  // std::cout << "shift_b = " << shift_b << std::endl;
  // std::cout << "split_value = " << split_value << std::endl;

  auto L_box = cur_mbr;
  auto R_box = cur_mbr;
  auto rx_prefix = x_prefix;
  auto ry_prefix = y_prefix;
  auto rz_prefix = z_prefix;
  if (split_axis == 0) {  //  cut x
    L_box.second.x = min(x_prefix + split_value - FT_EPS, L_box.second.x);
    R_box.first.x = max(x_prefix + split_value, R_box.first.x);
    rx_prefix += split_value;
  } else if (split_axis == 1) {  // cut y
    L_box.second.y = min(y_prefix + split_value - FT_EPS, L_box.second.y);
    R_box.first.y = max(y_prefix + split_value, R_box.first.y);
    ry_prefix += split_value;
  } else {  //  cut z
    L_box.second.z = min(z_prefix + split_value - FT_EPS, L_box.second.z);
    R_box.first.z = max(z_prefix + split_value, R_box.first.z);
    rz_prefix += split_value;
  }
  return make_tuple(L_box, R_box, rx_prefix, ry_prefix, rz_prefix);
}

auto get_delete_p(parlay::sequence<Point>& lhs, parlay::sequence<Point>& P,
                  size_t l, size_t r) {
  auto ret = parlay::sequence<Point>::uninitialized(lhs.size());
  auto pt_cmp = [&](auto& pt1, auto& pt2) {
    if (pt1.morton_id == pt2.morton_id) {
      if (pt1.id == pt2.id) {
        return 0;
      } else {
        if (pt1.id < pt2.id)
          return 1;
        else
          return 2;
      }
    } else {
      if (pt1.morton_id < pt2.morton_id)
        return 1;
      else
        return 2;
    }
  };
  size_t i = 0, j = l, cnt = 0;
  while (i < lhs.size() && j < r) {
    auto flag = pt_cmp(lhs[i], P[j]);
    if (!flag) {  // point should be deleted
      i++, j++;
    } else if (flag == 1) {  // first smaller, in A not in B
      parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
    } else {  //  second smaller, in B not in A
      j++;
    }
  }
  while (i < lhs.size()) {
    parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
  }
  ret.resize(cnt);
  return ret;
}

template <typename Func>
auto get_delete_p(parlay::sequence<Point>& lhs, parlay::sequence<Point>& P,
                  size_t l, size_t r, Func& f) {
  auto ret = parlay::sequence<Point>::uninitialized(lhs.size());
  auto pt_cmp = [&](auto& pt1, auto& pt2) {
    if (pt1.morton_id == pt2.morton_id) {
      if (pt1.id == pt2.id) {
        return 0;
      } else {
        if (pt1.id < pt2.id)
          return 1;
        else
          return 2;
      }
    } else {
      if (pt1.morton_id < pt2.morton_id)
        return 1;
      else
        return 2;
    }
  };
  size_t i = 0, j = l, cnt = 0;
  while (i < lhs.size() && j < r) {
    auto flag = pt_cmp(lhs[i], P[j]);
    if (!flag) {  // point should be deleted
      i++, j++;
    } else if (flag == 1) {  // first smaller, in A not in B
      if (!f(lhs[i])) {
        i++;
        continue;
      }
      parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
    } else {  //  second smaller, in B not in A
      j++;
    }
  }
  while (i < lhs.size()) {
    if (!f(lhs[i])) {
      i++;
      continue;
    }
    parlay::assign_uninitialized(ret[cnt++], lhs[i++]);
  }
  ret.resize(cnt);
  return ret;
}

// delete rhs from lhs, merge by morton_id and point id
template <typename T, typename P, typename DIFF>
auto merge_pts(P& a, T& b, DIFF& ret_diff, bool reverse = false) {
  size_t i = 0, j = 0;
  //  0: same point, 1: first is smaller, 2: second is smaller
  auto pt_cmp = [&](auto const& pt1, auto const& pt2) {
    if (pt1.morton_id == pt2.morton_id) {
      if (pt1.id == pt2.id) {
        return 0;
      } else {
        return pt1.id < pt2.id ? 1 : 2;
      }
    } else {
      return pt1.morton_id < pt2.morton_id ? 1 : 2;
    }
  };

  while (i < a.size() && j < b.size()) {
    auto flag = pt_cmp(a[i], b[j]);
    if (!flag) {  // same point
      i++, j++;
    } else if (flag == 1) {  // first smaller, in A not in B
      // parlay::assign_uninitialized(removed[cnt_remove++], a[i++]);
      ret_diff.remove_point(a[i], reverse);
      i++;
    } else {  //  second smaller, in B not in A
      // parlay::assign_uninitialized(added[cnt_add++], b[j++]);
      ret_diff.add_point(b[j], reverse);
      j++;
    }
  }
  while (i < a.size()) {
    // parlay::assign_uninitialized(removed[cnt_remove++], a[i++]);
    ret_diff.remove_point(a[i], reverse);
    i++;
  }
  while (j < b.size()) {
    // parlay::assign_uninitialized(added[cnt_add++], b[j++]);
    ret_diff.add_point(b[j], reverse);
    j++;
  }

  return;
}

// delete rhs from lhs, merge by morton_id and point id
template <typename T, typename P>
auto merge_pts(P& a, T& b) {
  size_t i = 0, j = 0;
  //  0: same point, 1: first is smaller, 2: second is smaller
  auto pt_cmp = [&](auto const& pt1, auto const& pt2) {
    if (pt1.morton_id == pt2.morton_id) {
      if (pt1.id == pt2.id) {
        return 0;
      } else {
        return pt1.id < pt2.id ? 1 : 2;
      }
    } else {
      return pt1.morton_id < pt2.morton_id ? 1 : 2;
    }
  };

  auto added = parlay::sequence<Point>::uninitialized(b.size());
  auto removed = parlay::sequence<Point>::uninitialized(a.size());

  size_t cnt_add = 0, cnt_remove = 0;

  while (i < a.size() && j < b.size()) {
    auto flag = pt_cmp(a[i], b[j]);
    if (!flag) {  // same point
      i++, j++;
    } else if (flag == 1) {  // first smaller, in A not in B
      parlay::assign_uninitialized(removed[cnt_remove++], a[i++]);
    } else {  //  second smaller, in B not in A
      parlay::assign_uninitialized(added[cnt_add++], b[j++]);
    }
  }
  while (i < a.size()) {
    parlay::assign_uninitialized(removed[cnt_remove++], a[i++]);
  }
  while (j < b.size()) {
    parlay::assign_uninitialized(added[cnt_add++], b[j++]);
  }
  removed.resize(cnt_remove);
  added.resize(cnt_add);

  return make_tuple(added, removed);
}

// delete rhs from lhs, merge by morton_id and point id, apply F
template <typename T, typename P, typename F, typename DIFF>
auto merge_pts(P& a, T& b, F& f, DIFF& ret_diff, bool reverse = false) {
  size_t i = 0, j = 0;
  //  0: same point, 1: first is smaller, 2: second is smaller
  auto pt_cmp = [&](auto const& pt1, auto const& pt2) {
    if (pt1.morton_id == pt2.morton_id) {
      if (pt1.id == pt2.id) {
        return 0;
      } else {
        return pt1.id < pt2.id ? 1 : 2;
      }
    } else {
      return pt1.morton_id < pt2.morton_id ? 1 : 2;
    }
  };

  while (i < a.size() && j < b.size()) {
    auto flag = pt_cmp(a[i], b[j]);

    if (!flag) {  // same point
      i++, j++;
    } else if (flag == 1) {  // first smaller, in A not in B
      if (f(a[i])) {
        ret_diff.remove_point(a[i], reverse);
        // parlay::assign_uninitialized(ret_diff.remove[ret_diff.remove_cnt++],
        // a[i]);
      }

      i++;
    } else {  //  second smaller, in B not in A
      if (f(b[j])) {
        ret_diff.add_point(b[j], reverse);
        // parlay::assign_uninitialized(ret_diff.add[ret_diff.add_cnt++], b[j]);
      }
      j++;
    }
  }
  while (i < a.size()) {
    if (f(a[i])) {
      ret_diff.remove_point(a[i], reverse);
      // parlay::assign_uninitialized(ret_diff.remove[ret_diff.remove_cnt++],
      // a[i]);
    }
    i++;
  }
  while (j < b.size()) {
    if (f(b[j])) {
      ret_diff.add_point(b[j], reverse);
      // parlay::assign_uninitialized(ret_diff.add[ret_diff.add_cnt++], b[j]);
    }
    j++;
  }

  return;
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
    // unsigned long a = P[i].x, b = P[i].y, c = P[i].z;
    // P[i].morton_id = libmorton::m3D_e_for<unsigned long long, unsigned
    // long>(a, b, c);
  });

  return parlay::sort(P, [&](auto lhs, auto rhs) {
    return lhs.morton_id < rhs.morton_id ||
           (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
  });
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

template <typename T>
auto diff_by_id(T& a, T& b) {
  auto id_cmp = [&](auto& lhs, auto& rhs) { return lhs.id < rhs.id; };

  auto sorted_a = parlay::sort(a, id_cmp);
  auto sorted_b = parlay::sort(b, id_cmp);

  auto added = T::uninitialized(b.size());
  auto removed = T::uninitialized(a.size());

  size_t cnt_add = 0, cnt_remove = 0;
  size_t i = 0, j = 0;

  while (i < a.size() && j < b.size()) {
    // cout << i << ", " << j << endl;
    if (sorted_a[i].id == sorted_b[j].id) {  // same point, check coordinate
      if (sorted_a[i] == sorted_b[j]) {      //  same coordinate, do nothing
        i++, j++;
      } else {
        parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
        parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
      }
    } else if (sorted_a[i].id <
               sorted_b[j].id) {  // first smaller, in A not in B
      parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
    } else {  //  second smaller, in B not in A
      parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
    }
  }
  while (i < sorted_a.size()) {
    parlay::assign_uninitialized(removed[cnt_remove++], sorted_a[i++]);
  }
  while (j < sorted_b.size()) {
    parlay::assign_uninitialized(added[cnt_add++], sorted_b[j++]);
  }
  removed.resize(cnt_remove);
  added.resize(cnt_add);

  return make_tuple(added, removed);
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
}  // namespace ZD3D
