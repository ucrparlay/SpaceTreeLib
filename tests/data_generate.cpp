#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "test_framework.h"

///**********************************START*********************************///

Coord const data_range = 1'000'000'000;

inline std::string toString(auto const& a) { return std::to_string(a); }

std::default_random_engine generator;
inline double getRealRandom(double const& a, double const& b) {
  std::uniform_real_distribution<double> distribution(a, b);
  return distribution(generator);
}

inline int getIntRandom(int const& a, int const& b) {
  std::uniform_int_distribution<int> distribution(a, b);
  return distribution(generator);
}

void print_to_file(std::string const& path, auto const& wp) {
  size_t const gen_num = wp.size();
  auto constexpr gen_dim =
      std::remove_reference_t<decltype(wp)>::value_type::GetDim();

  std::ofstream f;
  f.open(path);
  f << gen_num << " " << gen_dim << std::endl;

  for (size_t i = 0; i < gen_num; i++) {
    for (int j = 0; std::cmp_less(j, gen_dim); j++) {
      f << wp[i].pnt[j] << " ";  // TODO: format output by type
    }
    f << std::endl;
  }

  f.close();
}

template <typename Point>
class UniformGenerator {
 public:
  using Points = parlay::sequence<Point>;

  static Points WithinBox(size_t const gen_num) {
    Points wp(gen_num);

    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, data_range);

    parlay::random_generator gen(distrib(gen_mt));
    std::uniform_int_distribution<int> dis(
        0, data_range);  // WARN : assume using int

    // generate n random points in a cube
    parlay::parallel_for(0, gen_num, [&](size_t i) {
      auto r = gen[i];
      for (typename Point::DimsType j = 0; j < Point::GetDim(); j++) {
        wp[i].pnt[j] = dis(r);
      }
    });
    return wp;
  }

  static Points WithinSphere(size_t const gen_num, Point const& center,
                             typename Point::Coord const radius) {
    auto constexpr kDim = Point::GetDim();

    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, 65536);
    parlay::random_generator gen(
        distrib(gen_mt));  // PARA: thread safe random generator
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    return parlay::tabulate(gen_num, [&](size_t) -> Point {
      std::array<double, Point::GetDim()> coords;
      Point p;
      double r = static_cast<double>(radius) * std::pow(dis(gen), 1.0 / kDim);

      double sum_squares = 0.0;
      for (typename Point::DimsType j = 0; j < kDim; ++j) {
        coords[j] = dis(gen) * 2.0 - 1.0;
        sum_squares += coords[j] * coords[j];
      }

      double scale = r / std::sqrt(sum_squares);
      for (typename Point::DimsType j = 0; j < kDim; ++j) {
        coords[j] *= scale;
        p[j] = static_cast<typename Point::Coord>(coords[j]) + center[j];
      }

      assert(std::inner_product(coords.begin(), coords.end(), coords.begin(),
                                0) < radius * radius);

      return p;
    });
  }
};

template <typename Point, bool kSameDensity>
class VardenGenerator {
 public:
  using Points = parlay::sequence<Point>;
  using DimsType = Point::DimsType;

  constexpr static double GetRhoNoice() { return 1.0 / 10000; }

  constexpr static size_t GetCRest() { return 100; }

  static double GetRhoRestart(size_t n) {
    return 10.0 / (n * (1 - GetRhoNoice()));
  }

  static size_t GetNoicePtsNum(size_t gen_num) {
    return static_cast<size_t>(static_cast<double>(gen_num) * GetRhoNoice());
  }

  static size_t GetClusterPtsNum(size_t gen_num) {
    return gen_num - GetNoicePtsNum(gen_num);
  }

  size_t GetVincinity([[maybe_unused]] size_t const prev_restart) noexcept {
    if constexpr (kSameDensity) {
      return 100;
    } else {
      return 100 * ((prev_restart % 10) + 1);
    }
  }

  size_t GetShiftDistance([[maybe_unused]] size_t const prev_restart) noexcept {
    if constexpr (kSameDensity) {
      return 50 * Point::GetDim();
    } else {
      return GetVincinity(prev_restart) * Point::GetDim() / 2;
    }
  }

  Points GenerateCluster(size_t const gen_num, size_t const prev_restart) {
    using Spreader = std::pair<Point, DimsType>;
    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, data_range);

    parlay::random_generator gen(distrib(gen_mt));
    std::uniform_int_distribution<int> dis(0, data_range);

    auto generate_random_point = [&]() -> Point {
      Point p;
      for (auto& c : p.pnt) {
        c = dis(gen);
      }
      return p;
    };

    // ((a % b) + b) % b; negative modular
    auto shift_pts = parlay::tabulate(gen_num / GetCRest() + 1, [&](size_t i) {
      if (i) {
        int const dir = dis(gen) % Point::GetDim();
        Point p(static_cast<Point::Coord>(0));
        p[dir] = GetShiftDistance(prev_restart);
        return Spreader(p, dir);
      } else {
        return Spreader(generate_random_point(), 0);
      }
    });

    parlay::scan_inclusive_inplace(
        shift_pts.cut(1, shift_pts.size()),
        parlay::binary_op(
            [](Spreader const& a, Spreader const& b) {
              Point c(a.first);  // WARN: order here matters
              c[b.second] += b.first[b.second];
              return Spreader(c, b.second);
            },
            shift_pts[0]));

    // TODO: may need to adjust the range of coordinates
    auto pts =
        parlay::flatten(parlay::tabulate(shift_pts.size(), [&](size_t i) {
          size_t sz =
              i != shift_pts.size() - 1 ? GetCRest() : gen_num - i * GetCRest();
          return UniformGenerator<Point>::WithinSphere(
              sz, shift_pts[i].first, GetVincinity(prev_restart));
        }));

    assert(pts.size() == gen_num);
    return pts;
  }

  Points Apply(size_t const gen_num) {
    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, 2024);

    parlay::random_generator gen(
        distrib(gen_mt));  // PARA: thread safe random generator
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // NOTE: begin generate
    auto r_start_pts =
        parlay::tabulate(GetClusterPtsNum(gen_num), [&](size_t i) -> bool {
          auto r = gen[i];
          return i == 0 || dis(r) < GetRhoRestart(gen_num);
        });

    auto restart_idx = parlay::pack_index(r_start_pts);
    restart_idx.push_back(GetClusterPtsNum(gen_num));

    // for (auto i : restart_idx) {
    //   std::cout << i << " ";
    // }

    auto cluster_seq =
        parlay::flatten(parlay::tabulate(restart_idx.size() - 1, [&](size_t i) {
          return GenerateCluster(restart_idx[i + 1] - restart_idx[i], i);
        }));

    // std::cout << cluster_seq.size() << std::endl;

    return parlay::append(cluster_seq, UniformGenerator<Point>::WithinBox(
                                           GetNoicePtsNum(gen_num)));
  }

  size_t restarts_num;
  size_t r_vincinity;
  size_t r_shift;
};

template <typename Point>
auto generate_varden_points(size_t pts_num) {
  VardenGenerator<Point, false> varden;
  return varden.Apply(pts_num);
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>] [-s <kSummary>]");
  std::string path;
  int pts_num, pts_dim, file_num, varden;

  path = P.getOptionValue("-p");
  pts_num = P.getOptionIntValue("-n", 1'000'000);
  pts_dim = P.getOptionIntValue("-d", 2);
  file_num = P.getOptionIntValue("-f", 2);
  varden = P.getOptionIntValue("-t", 0);

  // NOTE: ../kdtree/ss_varden/
  path += (*path.rbegin() == '/' ? "" : "/") + toString(pts_num) + "_" +
          toString(pts_dim) + "/";
  std::filesystem::create_directory(path);

  auto generate = [&]<typename Point>(std::string const& new_path) {
    using Points = parlay::sequence<Point>;
    Points wp = varden ? VardenGenerator<Point, false>().Apply(pts_num)
                       : UniformGenerator<Point>::WithinBox(pts_num);
    print_to_file(new_path, wp);
  };

  for (int i = 0; i < file_num; i++) {
    std::string newpath = path + toString(i + 1) + ".in";
    std::cout << newpath << std::endl;
    if (pts_dim == 2) {
      generate.operator()<PointType<Coord, 2>>(newpath);
    } else if (pts_dim == 3) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else if (pts_dim == 5) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else if (pts_dim == 7) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else if (pts_dim == 9) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else if (pts_dim == 12) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else if (pts_dim == 16) {
      generate.operator()<PointType<Coord, 3>>(newpath);
    } else {
      throw std::runtime_error("Invalid dimension");
    }
  }
  return 0;
}
