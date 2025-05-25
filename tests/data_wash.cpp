#include "data_op.h"
#include "test_framework.h"

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
                "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
                "<_insertFile>] [-s <kSummary>]");
  std::string input_path, output_suffix, output_path;
  int pts_dim, multiply_offset;

  input_path = P.getOptionValue("-p");
  output_suffix = P.getOptionValue("-s", ".out");
  pts_dim = P.getOptionIntValue("-d", 2);
  multiply_offset = P.getOptionIntValue("-k", 1);

  std::cout << input_path << std::endl;
  std::cout << output_suffix << std::endl;

  auto generate = [&]<typename InputPoint>() {
    using OutputPoint = BasicPoint<int64_t, InputPoint::GetDim()>;
    using Laundy = DataLaundry<InputPoint, OutputPoint>;

    std::cout << "Input: " << input_path << std::endl;
    output_path = input_path.substr(0, input_path.rfind('.')) + output_suffix;
    std::cout << "Output: " << output_path << std::endl;

    parlay::sequence<InputPoint> wp;
    read_points<InputPoint>(input_path.c_str(), wp, 0);
    auto bb = pspt::BaseTree<InputPoint>::GetBox(parlay::make_slice(wp));
    std::cout << bb.first << " " << bb.second << std::endl;

    // PrintPoints(wp);
    auto new_wp = Laundy::RoundDown(wp, multiply_offset);
    // PrintPoints(new_wp);
    new_wp = Laundy::RemoveDuplicates(new_wp);
    // PrintPoints(new_wp);
    new_wp = Laundy::ShiftToFirstRegion(new_wp);
    // PrintPoints(new_wp);
    if (!Laundy::CheckCoordWithinRange(new_wp)) {
      throw std::runtime_error("ShiftToFirstRegion failed");
      exit(1);
    }

    std::cout << "Writing... " << std::endl;
    PrintToFile(output_path, new_wp);
  };

  if (pts_dim == 2) {
    generate.operator()<pspt::BasicPoint<double, 2>>();
  } else if (pts_dim == 3) {
    generate.operator()<pspt::BasicPoint<double, 3>>();
  } else if (pts_dim == 5) {
    generate.operator()<pspt::BasicPoint<double, 5>>();
  } else if (pts_dim == 7) {
    generate.operator()<pspt::BasicPoint<double, 7>>();
  } else if (pts_dim == 9) {
    generate.operator()<pspt::BasicPoint<double, 9>>();
  } else if (pts_dim == 12) {
    generate.operator()<pspt::BasicPoint<double, 12>>();
  } else if (pts_dim == 16) {
    generate.operator()<pspt::BasicPoint<double, 16>>();
  } else {
    throw std::runtime_error("Invalid dimension");
  }

  return 0;
}
