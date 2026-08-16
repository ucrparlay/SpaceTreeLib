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
#include "psi/base_tree.h"
#include "psi/dependence/basic_point.h"
#include "psi/dependence/comparator.h"
///**********************************START*********************************///

using axis_type = int64_t;
axis_type kValueUB = 1'000'000'000;
// axis_type const kValueUB = 1'000'000;

inline std::string toString(auto const &a)
{
	return std::to_string(a);
}

void print_points(auto const &wp)
{
	for (auto const &p : wp) {
		std::cout << p << std::endl;
	}
	std::cout << std::endl;
	return;
}

void print_to_file(std::string const &path, auto const &wp)
{
	size_t const gen_num = wp.size();
	auto constexpr gen_dim =
		std::remove_reference_t<decltype(wp)>::value_type::get_dim();

	auto header = toString(gen_num) + " " + toString(gen_dim);
	// NOTE: create a parlay::sequence<char> in order to call the
	// writeSeqToFile
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
		if (lineCount > 3)
			break;
	}
	if (lineCount <= 3) {
		throw std::runtime_error(
			"parallel written error, please retry." + path);
	}
	file.close();
	return;
}

template <typename InputPoint, typename output_point_type>
class data_laundry
{
public:
	using input_points_type = parlay::sequence<InputPoint>;
	using output_points_type = parlay::sequence<output_point_type>;

	static output_points_type round_down(auto const &wp,
					     int const multiply_offset)
	{
		return parlay::tabulate(
			wp.size(), [&](size_t i) -> output_point_type {
				output_point_type p;
				for (int j = 0;
				     j < output_point_type::get_dim(); j++) {
					p.pnt[j] = static_cast<
						typename output_point_type::
							coord_type>(
						std::floor(wp[i].pnt[j] *
							   multiply_offset));
				}
				p.aug = wp[i].aug;
				return p;
			});
	}

	static output_points_type remove_duplicates(auto const &wp)
	{
		using bp_type = typename output_point_type::bp_type;
		auto new_pts = parlay::unique(
			parlay::sort(
				wp,
				[&](auto const &a, auto const &b) {
					return static_cast<bp_type const &>(a) <
					       static_cast<bp_type const &>(b);
				}),
			[&](auto const &a, auto const &b) {
				return static_cast<bp_type const &>(a) ==
				       static_cast<bp_type const &>(b);
			});
		// return new_pts;
		return parlay::sort(new_pts, [&](auto const &a, auto const &b) {
			return a.aug.id < b.aug.id;
		});
	}

	static output_points_type shift_to_first_region(auto &wp)
	{
		auto bb = psi::base_tree<output_point_type>::get_box(
			parlay::make_slice(wp));
		return parlay::tabulate(wp.size(), [&](size_t i) {
			output_point_type p;
			for (int j = 0; j < output_point_type::get_dim(); j++) {
				p[j] = static_cast<
					typename output_point_type::coord_type>(
					wp[i][j] - bb.first[j]);
			}
			return p;
		});
	}

	static bool check_coord_within_range(auto &wp)
	{
		return parlay::all_of(wp, [&](auto const &p) {
			return std::ranges::all_of(p.pnt, [](auto const &c) {
				return c >= 0 && c <= INT32_MAX;
			});
		});
	}
};

template <typename Point>
class uniform_generator
{
public:
	using points_type = parlay::sequence<Point>;

	static points_type within_box(size_t const gen_num)
	{
		points_type wp(gen_num);

		std::random_device
			rd; // a seed source for the random number engine
		std::mt19937 gen_mt(
			rd()); // mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<int> distrib(1, kValueUB);

		parlay::random_generator gen(distrib(gen_mt));
		std::uniform_int_distribution<int> dis(
			0, kValueUB); // WARN : assume using int

		// generate n random points in a cube
		parlay::parallel_for(0, gen_num, [&](size_t i) {
			auto r = gen[i];
			for (typename Point::dims_type j = 0;
			     j < Point::get_dim(); j++) {
				wp[i].pnt[j] = dis(r);
			}
		});
		return wp;
	}

	static points_type
	within_sphere(size_t const gen_num, Point const &center,
		      typename Point::coord_type const radius)
	{
		auto constexpr num_dims = Point::get_dim();

		std::random_device
			rd; // a seed source for the random number engine
		std::mt19937 gen_mt(
			rd()); // mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<int> distrib(1, 65536);
		parlay::random_generator gen(
			distrib(gen_mt)); // PARA: thread safe random generator
		std::uniform_real_distribution<double> dis(0.0, 1.0);

		return parlay::tabulate(gen_num, [&](size_t) -> Point {
			std::array<double, Point::get_dim()> coords;
			Point p;
			double r = static_cast<double>(radius) *
				   std::pow(dis(gen), 1.0 / num_dims);

			double sum_squares = 0.0;
			for (typename Point::dims_type j = 0; j < num_dims;
			     ++j) {
				coords[j] = dis(gen) * 2.0 - 1.0;
				sum_squares += coords[j] * coords[j];
			}

			double scale = r / std::sqrt(sum_squares);
			for (typename Point::dims_type j = 0; j < num_dims;
			     ++j) {
				coords[j] *= scale;
				p[j] = static_cast<typename Point::coord_type>(
					       coords[j]) +
				       center[j];
			}

			assert(std::inner_product(coords.begin(), coords.end(),
						  coords.begin(),
						  0) < radius * radius);

			return p;
		});
	}
};

template <typename Point, bool SameDensity>
class varden_generator
{
public:
	using points_type = parlay::sequence<Point>;
	using dims_type = Point::dims_type;
	using num_type = psi::num_comparator<axis_type>;

	constexpr static double get_rho_noice()
	{
		return 1.0 / 10000;
	}

	// PARA: initial spreader avaliable counter
	constexpr static size_t spreader_size()
	{
		return 100;
	}

	static double restart_probability(size_t n)
	{
		return 10.0 / (n * (1 - get_rho_noice()));
	}

	static size_t get_noice_pts_num(size_t gen_num)
	{
		return static_cast<size_t>(static_cast<double>(gen_num) *
					   get_rho_noice());
	}

	static size_t get_cluster_pts_num(size_t gen_num)
	{
		return gen_num - get_noice_pts_num(gen_num);
	}

	size_t
	radius_by_vicinity([[maybe_unused]] size_t const prev_restart) noexcept
	{
		if constexpr (SameDensity) {
			return 100;
		} else {
			return 100 * ((prev_restart % 10) + 1);
		}
	}

	size_t
	shift_distance([[maybe_unused]] size_t const prev_restart) noexcept
	{
		if constexpr (SameDensity) {
			return 50 * Point::get_dim();
		} else {
			return radius_by_vicinity(prev_restart) *
			       Point::get_dim() / 2;
		}
	}

	points_type generate_cluster(size_t const gen_num,
				     size_t const cur_restart_idx)
	{
		using spreader_type = std::pair<Point, dims_type>;
		std::random_device
			rd; // a seed source for the random number engine
		std::mt19937 gen_mt(
			rd()); // mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<int> distrib(1, kValueUB);

		parlay::random_generator gen(distrib(gen_mt));
		std::uniform_int_distribution<int> dis(0, kValueUB);
		std::bernoulli_distribution bernoulli(0.5);

		auto generate_random_point = [&]() -> Point {
			Point p;
			for (auto &c : p.pnt) {
				c = dis(gen);
			}
			return p;
		};

		// ((a % b) + b) % b; negative modular
		auto shift_pts = parlay::tabulate(
			gen_num / spreader_size() + 1, [&](size_t i) {
				if (i) {
					int const dir =
						dis(gen) %
						Point::get_dim(); // pick a dir
					Point p(static_cast<Point::coord_type>(
						0));
					int sign = bernoulli(gen) ? 1 : -1;
					p[dir] =
						sign *
						shift_distance(
							cur_restart_idx); // move
									  // dis
									  // in
									  // this
									  // dir
					return spreader_type(p, dir);
				} else {
					return spreader_type(
						generate_random_point(), 0);
				}
			});

		parlay::scan_inclusive_inplace(
			shift_pts.cut(1, shift_pts.size()),
			parlay::binary_op(
				[](spreader_type const &a,
				   spreader_type const &b) {
					Point c(a.first); // WARN: order here
							  // matters
					c[b.second] += b.first[b.second];
					return spreader_type(c, b.second);
				},
				shift_pts[0]));

		// TODO: may need to adjust the range of coordinates
		auto pts = parlay::flatten(
			parlay::tabulate(shift_pts.size(), [&](size_t i) {
				size_t sz =
					i != shift_pts.size() - 1
						? spreader_size()
						: gen_num - i * spreader_size();
				return uniform_generator<Point>::within_sphere(
					sz, shift_pts[i].first,
					radius_by_vicinity(cur_restart_idx));
			}));

		assert(pts.size() == gen_num);
		return pts;
	}

	points_type apply(size_t const gen_num)
	{
		std::random_device
			rd; // a seed source for the random number engine
		std::mt19937 gen_mt(
			rd()); // mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<int> distrib(1, 2024);

		parlay::random_generator gen(
			distrib(gen_mt)); // PARA: thread safe random generator
		std::uniform_real_distribution<double> dis(0.0, 1.0);

		// NOTE: begin generate
		// pick the restart points
		auto r_start_pts = parlay::tabulate(
			get_cluster_pts_num(gen_num), [&](size_t i) -> bool {
				auto r = gen[i];
				return i == 0 ||
				       dis(r) < restart_probability(gen_num);
			});

		// pack the index of the restart points
		auto restart_idx = parlay::pack_index(r_start_pts);
		restart_idx.push_back(get_cluster_pts_num(gen_num));

		// generate the points for each cluster and pack them
		auto cluster_seq = parlay::flatten(
			parlay::tabulate(restart_idx.size() - 1, [&](size_t i) {
				return generate_cluster(
					restart_idx[i + 1] - restart_idx[i], i);
			}));

		// move the points to the first region of the space
		auto bb = psi::base_tree<Point>::get_box(
			parlay::make_slice(cluster_seq));
		for (dims_type j = 0; j < Point::get_dim(); j++) {
			bb.first[j] = num_type::min(bb.first[j],
						    static_cast<axis_type>(0));
			bb.first[j] = num_type::abs(bb.first[j]);
		}
		parlay::parallel_for(0, cluster_seq.size(), [&](size_t i) {
			for (dims_type j = 0; j < Point::get_dim(); j++) {
				cluster_seq[i][j] += bb.first[j];
				if (cluster_seq[i][j] < 0) {
					std::cout << "cluster_seq[i][j]: "
						  << cluster_seq[i][j]
						  << std::endl;
				}
			}
		});

		// append the noice points
		return parlay::append(cluster_seq,
				      uniform_generator<Point>::within_box(
					      get_noice_pts_num(gen_num)));
	}
};

template <typename Point>
auto generate_varden_points(size_t pts_num)
{
	varden_generator<Point, false> varden;
	return varden.apply(pts_num);
}
