#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>

#include "common/IO.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pspt/base_tree.h"
#include "pspt/dependence/basic_point.h"
#include "pspt/dependence/comparator.h"
///**********************************START*********************************///

using Axis = int64_t;
Axis const kValueUB = 1'000'000'000;
// Coord const kValueUB = 1'000'000;

inline std::string toString(auto const& a) { return std::to_string(a); }

void PrintToFile(std::string const& path, auto const& wp) {
  size_t const gen_num = wp.size();
  auto constexpr gen_dim =
      std::remove_reference_t<decltype(wp)>::value_type::GetDim();

  auto header = toString(gen_num) + " " + toString(gen_dim);
  // NOTE: create a parlay::sequence<char> in order to call the writeSeqToFile
  auto wp_str = parlay::tabulate(wp.size(), [&](size_t i) {
    std::string str;
    for (int j = 0; std::cmp_less(j, gen_dim); j++) {
      str += std::to_string(wp[i].pnt[j]) + " ";
    }
    parlay::sequence<char> cstr(str.begin(), str.end());
    return cstr;
  });

  // NOTE: write to file in parallel
  writeSeqToFile(header, wp_str, path.c_str());

  // NOTE: sometimes parallel written cannot be finished correcly
  // naive check for it
  std::ifstream file(path);
  int lineCount = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++lineCount;
    if (lineCount > 3) break;
  }
  if (lineCount <= 3) {
    throw std::runtime_error("parallel written error, please retry." + path);
  }
  file.close();
  return;
}

template <typename Point>
class UniformGenerator {
 public:
  using Points = parlay::sequence<Point>;

  static Points WithinBox(size_t const gen_num) {
    Points wp(gen_num);

    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, kValueUB);

    parlay::random_generator gen(distrib(gen_mt));
    std::uniform_int_distribution<int> dis(
        0, kValueUB);  // WARN : assume using int

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
  using Num = pspt::Num_Comparator<Axis>;

  constexpr static double GetRhoNoice() { return 1.0 / 10000; }

  // PARA: initial spreader avaliable counter
  constexpr static size_t SpreaderSize() { return 100; }

  static double RestartProbability(size_t n) {
    return 10.0 / (n * (1 - GetRhoNoice()));
  }

  static size_t GetNoicePtsNum(size_t gen_num) {
    return static_cast<size_t>(static_cast<double>(gen_num) * GetRhoNoice());
  }

  static size_t GetClusterPtsNum(size_t gen_num) {
    return gen_num - GetNoicePtsNum(gen_num);
  }

  size_t RadisuByVincinity(
      [[maybe_unused]] size_t const prev_restart) noexcept {
    if constexpr (kSameDensity) {
      return 100;
    } else {
      return 100 * ((prev_restart % 10) + 1);
    }
  }

  size_t ShiftDistance([[maybe_unused]] size_t const prev_restart) noexcept {
    if constexpr (kSameDensity) {
      return 50 * Point::GetDim();
    } else {
      return RadisuByVincinity(prev_restart) * Point::GetDim() / 2;
    }
  }

  Points GenerateCluster(size_t const gen_num, size_t const cur_restart_idx) {
    using Spreader = std::pair<Point, DimsType>;
    std::random_device rd;      // a seed source for the random number engine
    std::mt19937 gen_mt(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(1, kValueUB);

    parlay::random_generator gen(distrib(gen_mt));
    std::uniform_int_distribution<int> dis(0, kValueUB);
    std::bernoulli_distribution bernoulli(0.5);

    auto generate_random_point = [&]() -> Point {
      Point p;
      for (auto& c : p.pnt) {
        c = dis(gen);
      }
      return p;
    };

    // ((a % b) + b) % b; negative modular
    auto shift_pts =
        parlay::tabulate(gen_num / SpreaderSize() + 1, [&](size_t i) {
          if (i) {
            int const dir = dis(gen) % Point::GetDim();  // pick a dir
            Point p(static_cast<Point::Coord>(0));
            int sign = bernoulli(gen) ? 1 : -1;
            p[dir] =
                sign * ShiftDistance(cur_restart_idx);  // move dis in this dir
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
          size_t sz = i != shift_pts.size() - 1 ? SpreaderSize()
                                                : gen_num - i * SpreaderSize();
          return UniformGenerator<Point>::WithinSphere(
              sz, shift_pts[i].first, RadisuByVincinity(cur_restart_idx));
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
    // pick the restart points
    auto r_start_pts =
        parlay::tabulate(GetClusterPtsNum(gen_num), [&](size_t i) -> bool {
          auto r = gen[i];
          return i == 0 || dis(r) < RestartProbability(gen_num);
        });

    // pack the index of the restart points
    auto restart_idx = parlay::pack_index(r_start_pts);
    restart_idx.push_back(GetClusterPtsNum(gen_num));

    // generate the points for each cluster and pack them
    auto cluster_seq =
        parlay::flatten(parlay::tabulate(restart_idx.size() - 1, [&](size_t i) {
          return GenerateCluster(restart_idx[i + 1] - restart_idx[i], i);
        }));

    // move the points to the first region of the space
    auto bb = pspt::BaseTree<Point>::GetBox(parlay::make_slice(cluster_seq));
    for (DimsType j = 0; j < Point::GetDim(); j++) {
      bb.first[j] = Num::Min(bb.first[j], static_cast<Axis>(0));
      bb.first[j] = Num::Abs(bb.first[j]);
    }
    parlay::parallel_for(0, cluster_seq.size(), [&](size_t i) {
      for (DimsType j = 0; j < Point::GetDim(); j++) {
        cluster_seq[i][j] += bb.first[j];
        if (cluster_seq[i][j] < 0) {
          std::cout << "cluster_seq[i][j]: " << cluster_seq[i][j] << std::endl;
        }
      }
    });

    // append the noice points
    return parlay::append(cluster_seq, UniformGenerator<Point>::WithinBox(
                                           GetNoicePtsNum(gen_num)));
  }
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
  std::cout << "pts_num: " << pts_num << " pts_dim: " << pts_dim
            << " file_num: " << file_num << " varden: " << varden << std::endl;

  // NOTE: ../kdtree/
  path += std::string(*path.rbegin() == '/' ? "" : "/") +
          std::string(varden ? "ss_varden_bigint/" : "uniform_bigint/") +
          toString(pts_num) + "_" + toString(pts_dim) + "/";
  // + "/0_" + toString(static_cast<size_t>(kValueUB)) + "/";
  std::filesystem::create_directory(path);

  auto generate = [&]<typename Point>(std::string const& new_path) {
    std::cout << "Generating... " << std::endl;
    auto wp = varden ? VardenGenerator<Point, false>().Apply(pts_num)
                     : UniformGenerator<Point>::WithinBox(pts_num);
    std::cout << "Writing... " << std::endl;
    PrintToFile(new_path, wp);
  };

  for (int i = 0; i < file_num; i++) {
    std::string newpath = path + toString(i + 1) + ".in";
    std::cout << newpath << std::endl;
    if (pts_dim == 2) {
      generate.operator()<pspt::BasicPoint<Axis, 2>>(newpath);
    } else if (pts_dim == 3) {
      generate.operator()<pspt::BasicPoint<Axis, 3>>(newpath);
    } else if (pts_dim == 5) {
      generate.operator()<pspt::BasicPoint<Axis, 5>>(newpath);
    } else if (pts_dim == 7) {
      generate.operator()<pspt::BasicPoint<Axis, 7>>(newpath);
    } else if (pts_dim == 9) {
      generate.operator()<pspt::BasicPoint<Axis, 9>>(newpath);
    } else if (pts_dim == 12) {
      generate.operator()<pspt::BasicPoint<Axis, 12>>(newpath);
    } else if (pts_dim == 16) {
      generate.operator()<pspt::BasicPoint<Axis, 16>>(newpath);
    } else {
      throw std::runtime_error("Invalid dimension");
    }
  }
  return 0;
}
