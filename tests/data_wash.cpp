#include "data_op.h"
#include "test_framework.h"

template <typename InputPoint>
void round_data_and_shift(commandLine &P)
{
	std::string input_path, output_suffix, output_path;
	int pts_dim, multiply_offset;

	input_path = P.getOptionValue("-p");
	output_suffix = P.getOptionValue("-output_suffix", ".out");
	pts_dim = P.getOptionIntValue("-d", 2);
	multiply_offset = P.getOptionIntValue("-multiply_offset", 1);

	std::cout << input_path << std::endl;
	std::cout << output_suffix << std::endl;
	using output_point_type =
		aug_point<int64_t, InputPoint::get_dim(), wrapper::aug_id>;
	using laundy_type = data_laundry<InputPoint, output_point_type>;

	std::cout << "Input: " << input_path << std::endl;
	output_path =
		input_path.substr(0, input_path.rfind('.')) + output_suffix;
	std::cout << "Output: " << output_path << std::endl;

	parlay::sequence<InputPoint> wp;
	read_points<InputPoint>(input_path.c_str(), wp, 0);
	auto bb = psi::base_tree<InputPoint>::get_box(parlay::make_slice(wp));
	std::cout << bb.first << " " << bb.second << std::endl;

	// print_points(wp);
	auto new_wp = laundy_type::round_down(wp, multiply_offset);
	// print_points(new_wp);
	new_wp = laundy_type::remove_duplicates(new_wp);
	// print_points(new_wp);
	new_wp = laundy_type::shift_to_first_region(new_wp);
	// print_points(new_wp);
	if (!laundy_type::check_coord_within_range(new_wp)) {
		throw std::runtime_error("shift_to_first_region failed");
		exit(1);
	}

	auto new_bb = psi::base_tree<output_point_type>::get_box(
		parlay::make_slice(new_wp));
	std::cout << new_bb.first << " " << new_bb.second << std::endl;

	std::cout << "Writing... " << std::endl;
	print_to_file(output_path, new_wp);
}

template <typename InputPoint>
void sort_by_coord(commandLine &P)
{
	std::string input_path, output_suffix, output_path;
	int pts_dim, sort_dim;

	input_path = P.getOptionValue("-p");
	output_suffix = P.getOptionValue("-output_suffix", ".out");
	pts_dim = P.getOptionIntValue("-d", 2);
	sort_dim = P.getOptionIntValue("-sort_dim", 0);

	using output_point_type =
		aug_point<int64_t, InputPoint::get_dim(), wrapper::aug_id>;
	using laundy_type = data_laundry<InputPoint, output_point_type>;

	std::cout << "Input: " << input_path << std::endl;
	output_path =
		input_path.substr(0, input_path.rfind('.')) + output_suffix;
	std::cout << "Output: " << output_path << std::endl;

	parlay::sequence<InputPoint> wp;
	read_points<InputPoint>(input_path.c_str(), wp, 0);

	auto new_wp = parlay::sort(wp, [&](auto const &a, auto const &b) {
		return a.pnt[sort_dim] < b.pnt[sort_dim];
	});

	std::cout << "Writing... " << std::endl;
	print_to_file(output_path, new_wp);
}

int main(int argc, char *argv[])
{
	commandLine P(
		argc, argv,
		"[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
		"<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
		"<_insertFile>] [-s <kSummary>] [-u <usage>]");
	int coord_type = P.getOptionIntValue("-coord_type", 0);
	int pts_dim = P.getOptionIntValue("-d", 2);
	int usage = P.getOptionIntValue("-usage", 0);

	auto apply = [&]<typename InputPoint>() {
		if (usage & (1 << 0)) {
			round_data_and_shift<InputPoint>(P);
		} else if (usage & (1 << 1)) {
			sort_by_coord<InputPoint>(P);
		}
	};

	auto generate_with_coord_type = [&]<typename CoordType>() {
		if (pts_dim == 2) {
			apply.template operator()<psi::aug_point<
				CoordType, 2, wrapper::aug_id>>();
		} else if (pts_dim == 3) {
			apply.template operator()<psi::aug_point<
				CoordType, 3, wrapper::aug_id>>();
		} else if (pts_dim == 5) {
			apply.template operator()<psi::aug_point<
				CoordType, 5, wrapper::aug_id>>();
		} else if (pts_dim == 7) {
			apply.template operator()<psi::aug_point<
				CoordType, 7, wrapper::aug_id>>();
		} else if (pts_dim == 9) {
			apply.template operator()<psi::aug_point<
				CoordType, 9, wrapper::aug_id>>();
		} else if (pts_dim == 12) {
			apply.template operator()<psi::aug_point<
				CoordType, 12, wrapper::aug_id>>();
		} else if (pts_dim == 16) {
			apply.template operator()<psi::aug_point<
				CoordType, 16, wrapper::aug_id>>();
		} else {
			throw std::runtime_error("Invalid dimension");
		}
	};

	// Decide coordinate type based on coord_type parameter
	if (coord_type == 0) {
		generate_with_coord_type.template operator()<long>();
	} else if (coord_type == 1) {
		generate_with_coord_type.template operator()<double>();
	} else {
		throw std::runtime_error("Invalid coordinate type: " +
					 std::to_string(coord_type));
	}

	return 0;
}
