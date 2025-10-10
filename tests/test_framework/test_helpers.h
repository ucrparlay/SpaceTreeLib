#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>

#include "../common/geometryIO.h"
#include "baselines/zdtree/zdtree.hpp"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/slice.h"
#include "psi/dependence/basic_point.h"
#include "psi/dependence/type_trait.h"
#include "test_config.h"

// Helper function for generating rectangles
template <typename Point, typename Tree, bool SavePoint>
size_t recurse_box(parlay::slice<Point*, Point*> In, auto& box_seq, int DIM,
                   std::pair<size_t, size_t> range, int& idx, int rec_num,
                   int type) {
  using Geo = Tree::Geo;
  size_t n = In.size();
  if (idx >= rec_num || n < range.first || n == 0) return 0;

  size_t mx = 0;
  bool goon = false;
  if (n <= range.second) {
    if constexpr (SavePoint) {
      box_seq[idx++] = std::make_pair(Geo::GetBox(In), parlay::to_sequence(In));
    } else {
      box_seq[idx++] = std::make_pair(Geo::GetBox(In), In.size());
    }

    if (parlay::all_of(
            In, [&](Point const& p) { return p.SameDimension(In[0]); })) {
      // NOTE: handle the case that all Points are the same which is
      // undivideable
      return In.size();
    } else {
      goon = true;
      mx = n;
    }
  }

  int dim = get_random_index(0, DIM, rand());
  size_t pos = get_random_index(0, n, rand());
  parlay::sequence<bool> flag(n, 0);
  parlay::parallel_for(0, n, [&](size_t i) {
    if (psi::Num_Comparator<Coord>::Gt(In[i][dim], In[pos][dim]))
      flag[i] = 1;
    else
      flag[i] = 0;
  });
  auto [Out, m] = parlay::internal::split_two(In, flag);

  assert(Out.size() == n);
  size_t l = 0, r = 0;
  l = recurse_box<Point, Tree, SavePoint>(Out.cut(0, m), box_seq, DIM, range,
                                          idx, rec_num, type);
  r = recurse_box<Point, Tree, SavePoint>(Out.cut(m, n), box_seq, DIM, range,
                                          idx, rec_num, type);

  if (goon) {
    return mx;
  } else {
    return std::max(l, r);
  }
}

template <typename Point, typename Tree, bool SavePoint, bool FixSize = false>
auto gen_rectangles(int rec_num, int const type,
                    parlay::sequence<Point> const& WP, int DIM) {
  using Points = typename Tree::Points;
  using Box = typename Tree::Box;
  using BoxSeq = std::conditional_t<SavePoint == false,
                                    parlay::sequence<std::pair<Box, size_t>>,
                                    parlay::sequence<std::pair<Box, Points>>>;

  size_t n = WP.size();
  std::pair<size_t, size_t> range;
  if constexpr (FixSize) {
    if (type == 0) {  //* small bracket
      range.first = 1;
      range.second = 100;
    } else if (type == 1) {  //* medium bracket
      range.first = 100;
      range.second = 10000;
    } else if (type == 2) {  //* large bracket
      range.first = 10000;
      if (n > 1'000'000) {
        range.second = 1'000'000;
      }  // ensure we can generate 50k rect.
      else {
        range.second = n - 1;
      }
    }
  } else {
    if (type == 0) {  //* small bracket
      range.first = 1;
      range.second = static_cast<size_t>(std::sqrt(std::sqrt(1.0 * n)));
    } else if (type == 1) {  //* medium bracket
      range.second = static_cast<size_t>(std::sqrt(std::sqrt(1.0 * n)));
      range.second = static_cast<size_t>(std::sqrt(1.0 * n));
    } else if (type == 2) {  //* large bracket
      range.second = static_cast<size_t>(std::sqrt(1.0 * n));
      if (n > 1'000'000) {
        range.second = 1'000'000;
      }  // ensure we can generate 50k rect.
      else {
        range.second = n - 1;
      }
    }
  }
  BoxSeq box_seq(rec_num);
  int cnt = 0;
  Points wp(n);

  srand(10);

  size_t max_size = 0;
  while (cnt < rec_num) {
    parlay::copy(WP, wp);
    auto r = recurse_box<Point, Tree, SavePoint>(
        parlay::make_slice(wp), box_seq, DIM, range, cnt, rec_num, type);
    max_size = std::max(max_size, r);
  }
  return std::make_pair(box_seq, max_size);
}

//*---------- generate Points within a 0-box_size --------------------
template <typename Point>
void generate_random_points(parlay::sequence<Point>& wp, Coord _box_size,
                            long n, int Dim) {
  Coord box_size = _box_size;

  std::random_device rd;      // a seed source for the random number engine
  std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distrib(1, box_size);

  parlay::random_generator gen(distrib(gen_mt));
  std::uniform_int_distribution<int> dis(0, box_size);

  wp.resize(n);
  // generate n random Points in a cube
  parlay::parallel_for(
      0, n,
      [&](long i) {
        auto r = gen[i];
        for (int j = 0; j < Dim; j++) {
          wp[i][j] = dis(r);
        }
      },
      1000);
  return;
}

template <typename Point>
std::pair<size_t, int> read_points(char const* iFile,
                                   parlay::sequence<Point>& wp,
                                   [[maybe_unused]] int id_offset = 0) {
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  static Coords a_sample_point;
  parlay::sequence<char> S = readStringFromFile(iFile);
  parlay::sequence<char*> W = stringToWords(S);
  size_t N = std::stoul(W[0], nullptr, 10);
  int Dim = atoi(W[1]);
  assert(N > 0 && Dim >= 1);

  auto pts = W.cut(2, W.size());
  assert(pts.size() % Dim == 0);
  size_t n = pts.size() / Dim;
  auto a = parlay::tabulate(Dim * n, [&](size_t i) -> Coord {
    if constexpr (std::is_integral_v<Coord>)
      return std::stol(pts[i]);
    else if (std::is_floating_point_v<Coord>)
      return std::stod(pts[i]);
  });
  wp.resize(N);
  parlay::parallel_for(0, n, [&](size_t i) {
    for (int j = 0; j < Dim; j++) {
      wp[i][j] = a[i * Dim + j];
      if constexpr (std::is_same_v<
                        Point, psi::BasicPoint<Coord, a_sample_point.size()>>) {
        ;
      } else if constexpr (std::is_same_v<Point, typename ZD::geobase::Point>) {
        wp[i].id = i + id_offset;
      } else {
        wp[i].aug.id = i + id_offset;
      }
    }
  });
  return std::make_pair(N, Dim);
}

template <typename T>
class counter_iterator {
 private:
  struct accept_any {
    template <typename U>
    accept_any& operator=(U const&) {
      return *this;
    }
  };

 public:
  typedef std::output_iterator_tag iterator_category;

  counter_iterator(T& counter) : counter(counter) {}
  counter_iterator& operator=(counter_iterator const& other) {
    if (this != &other) {  // Check for self-assignment
      counter = other.counter;
    }
    return *this;
  }
  counter_iterator(counter_iterator const& other) : counter(other.counter) {}

  bool operator==(counter_iterator const& rhs) const {
    return counter == rhs.counter;
  }
  bool operator!=(counter_iterator const& rhs) const {
    return counter != rhs.counter;
  }

  accept_any operator*() const {
    ++counter.get();
    return {};
  }

  counter_iterator& operator++() {  // ++a
    return *this;
  }
  counter_iterator operator++(int) {  // a++
    return *this;
  }

 protected:
  std::reference_wrapper<T> counter;
};

struct StepUpdateLogger {
  int id;
  double t;

  friend std::ostream& operator<<(std::ostream& os,
                                  StepUpdateLogger const& log) {
    os << "(" << log.id << ", " << log.t << ")";
    return os;
  }
};
