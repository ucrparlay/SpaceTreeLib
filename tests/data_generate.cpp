#include <filesystem>
#include <iostream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "test_framework.h"

///**********************************START*********************************///
std::string path;
int pts_num, pts_dim, file_num, varden;

Coord const data_range = 1e9;

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
  std::ofstream f;
  f.open(path);
  f << pts_num << " " << pts_dim << std::endl;
  for (size_t i = 0; std::cmp_less(i, pts_num); i++) {
    for (int j = 0; j < pts_dim; j++) {
      f << wp[i].pnt[j] << " ";
    }
    f << std::endl << std::flush;
  }
  f.close();
}

template <typename Point, bool kSameDensity>
struct VardenParam {
  using Points = parlay::sequence<Point>;

  VardenParam()
    requires(kSameDensity)
      : i(0), r_vincinity(100), r_shift(50 * Point::GetDim()) {}

  VardenParam()
    requires(!kSameDensity)
      : i(0),
        r_vincinity(100 * ((i % 10) + 1)),
        r_shift(r_vincinity * Point::GetDim() / 2) {}

  constexpr static double GetRhoNoice() { return 1.0 / 10000; }
  constexpr static size_t GetRest() { return 100; }
  static double GetRhoRestart(size_t n) {
    return 10.0 / (n * (1 - GetRhoNoice()));
  }
  void UpdateVincinity() noexcept { r_vincinity = 100 * ((i % 10) + 1); }

  Points Apply() {
    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, 2024);

    parlay::random_generator gen(distrib(gen_mt));
    std::uniform_real_distribution<double> dis(0, 1);

    auto r_start_pts = parlay::tabulate(pts_num, [&](size_t i) -> bool {
      auto r = gen[i];
      return i == 0 || dis(r) < GetRhoRestart(pts_num);
    });
    // auto reduce_pts =
    //     parlay::scan(r_start_pts, parlay::binary_op(std::plus<size_t>(), 0));
    auto restart_idx = parlay::pack_index(r_start_pts);
    restart_idx.push_back(pts_num);
    return parlay::flatten(
        parlay::tabulate(restart_idx.size() - 1, [&](size_t i) {
          auto l = restart_idx[i];
          auto r = restart_idx[i + 1];
          Points cluster;
          return cluster;
        }));
  }

  size_t i;
  size_t r_vincinity;
  size_t r_shift;
};

template <typename Point>
auto generate_uniform_points() {
  using Points = parlay::sequence<Point>;
  Points wp(pts_num);

  std::random_device rd;      // a seed source for the random number engine
  std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distrib(1, data_range);

  parlay::random_generator gen(distrib(gen_mt));
  std::uniform_int_distribution<int> dis(0, data_range);

  // generate n random points in a cube
  parlay::parallel_for(
      0, pts_num,
      [&](size_t i) {
        auto r = gen[i];
        for (int j = 0; j < pts_dim; j++) {
          wp[i].pnt[j] = dis(r);
        }
      },
      1000);
  return wp;
}

template <typename Point>
auto generate_varden_points() {
  VardenParam<Point, false> varden;
  return varden.Apply();
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>] [-s <kSummary>]");
  path = P.getOptionValue("-p");
  pts_num = P.getOptionIntValue("-n", 1'000'000);
  pts_dim = P.getOptionIntValue("-d", 2);
  file_num = P.getOptionIntValue("-f", 2);
  varden = P.getOptionIntValue("-t", 0);

  path += "/" + toString(pts_num) + "_" + toString(pts_dim) + "/";
  std::filesystem::create_directory(path);

  auto generate = [&]<typename Point>(std::string const& new_path) {
    using Points = parlay::sequence<Point>;
    Points wp = varden ? generate_varden_points<Point>()
                       : generate_uniform_points<Point>();
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
