#include "data_op.h"

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
  file_num = P.getOptionIntValue("-file_num", 2);
  varden = P.getOptionIntValue("-varden", 0);
  std::cout << "pts_num: " << pts_num << " pts_dim: " << pts_dim
            << " file_num: " << file_num << " varden: " << varden << std::endl;

  // NOTE: ../kdtree/
  path += std::string(*path.rbegin() == '/' ? "" : "/") +
          // std::string(varden ? "ss_varden_bigint/" : "uniform_bigint/") +
          std::string(varden ? "ss_varden/" : "uniform/") + toString(pts_num) +
          "_" + toString(pts_dim) + "/";
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
      generate.operator()<psi::BasicPoint<Axis, 2>>(newpath);
    } else if (pts_dim == 3) {
      generate.operator()<psi::BasicPoint<Axis, 3>>(newpath);
    } else if (pts_dim == 5) {
      generate.operator()<psi::BasicPoint<Axis, 5>>(newpath);
    } else if (pts_dim == 7) {
      generate.operator()<psi::BasicPoint<Axis, 7>>(newpath);
    } else if (pts_dim == 9) {
      generate.operator()<psi::BasicPoint<Axis, 9>>(newpath);
    } else if (pts_dim == 12) {
      generate.operator()<psi::BasicPoint<Axis, 12>>(newpath);
    } else if (pts_dim == 16) {
      generate.operator()<psi::BasicPoint<Axis, 16>>(newpath);
    } else {
      throw std::runtime_error("Invalid dimension");
    }
  }
  return 0;
}
