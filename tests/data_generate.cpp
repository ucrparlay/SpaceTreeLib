#include <filesystem>
#include <iostream>
#include <string>

#include "test_framework.h"

constexpr int MAX_DIM = 16;

///**********************************START*********************************///
std::string path;
int pts_num, pts_dim, file_num, distribution;

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

template <typename Point>
void generate_uniform_points(auto& wp) {
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
  return;
}

template <typename Point>
void generate_varden_points(auto& wp) {
  assert(wp.size());
  return;
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
  distribution = P.getOptionIntValue("-t", 0);

  path += "/" + toString(pts_num) + "_" + toString(pts_dim) + "/";
  std::filesystem::create_directory(path);

  auto generate = [&]<typename Point>(std::string const& new_path) {
    using Points = parlay::sequence<Point>;
    Points wp(pts_num);
    if (distribution == 0) {
      generate_uniform_points<Point>(wp);
    } else {
      generate_varden_points<Point>(wp);
    }
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
