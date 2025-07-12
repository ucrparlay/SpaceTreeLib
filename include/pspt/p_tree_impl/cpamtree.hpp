#pragma once 

#include <bits/stdc++.h>

#include <cpam/cpam.h>
// #include <pam/pam.h>
#include <parlay/primitives.h>
#include "geobase.h"
#include <parlay/internal/get_time.h>
#include "pam/utils.h"

namespace CPAMTree{
	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	using key_type = pair<unsigned long long, unsigned long long>;	// morton_id, id
	using val_type = Point;
	using aug_type = pair<Bounding_Box, size_t>;

	//	CPAM entry
	struct entry {
		using key_t = key_type;
		using val_t = val_type;
		using aug_t = aug_type;

		static inline bool comp(key_t a, key_t b) { return a < b; }
		static aug_t get_empty() { return make_pair(Bounding_Box{Point(1e60, 1e60), Point(-1, -1)}, 0); }
		static aug_t from_entry(key_t k, val_t v) { return make_pair(Bounding_Box(v, v), 1); }
		static aug_t combine(aug_t a, aug_t b) { return make_pair(merge_mbr(a.first, b.first), a.second + b.second); }
	};

	using zmap = cpam::aug_map<entry, 32>;
	using par = std::tuple<entry::key_t, entry::val_t>;

	template<typename T>
	auto knn(T &tree, geobase::Point &query_point, size_t &k){
		auto f = [&](auto cur_pt){ return point_point_sqrdis(cur_pt, query_point); };

		auto f2 = [&](auto cur_mbr){ return point_mbr_sqrdis(query_point, cur_mbr); };

		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		zmap::knn_filter(tree, f, f2, k, nn_res);
		return nn_res.top().second;
	}


	template<typename M>
	auto map_diff(M &lhs, M &rhs){
		auto add = zmap::values(zmap::map_difference(rhs, lhs));
		auto remove = zmap::values(zmap::map_difference(lhs, rhs));
		return make_tuple(add, remove);
	}

	template<class PT>
	auto map_init(PT &P, bool use_hilbert = false){
		size_t n = P.size();
		
		parlay::internal::timer t("SFC time", true);
		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		// t.stop();
		// t.total();
		// auto SFC_time = t.total_time();
		auto h_values = parlay::sequence<unsigned long long>::uninitialized(n);
		auto h_values2 =parlay::sequence<unsigned long long>::uninitialized(n); 
		parlay::internal::timer t1("Create entry time", true);
		parlay::sequence<par> entries(n);
		parlay::parallel_for(0, n, [&](int i){
			entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// entries[i] = {P[i]->id, P[i]};
			h_values[i] = P[i].morton_id;
			h_values2[i] = P[i].morton_id;
		});
		// t1.total();
		
		parlay::internal::timer t_integer_sort("integer_sort", true);
		auto ret = parlay::integer_sort(h_values);
		// t_integer_sort.total();
		auto less = [&](int a, int b){
			return a < b;
		};
		parlay::internal::timer t_sample_sort("sample_sort", true);
		auto ret2 = parlay::internal::sample_sort(parlay::make_slice(h_values2.begin(), h_values2.end()), less);
		// t_sample_sort.total();

		parlay::internal::timer t2("CPAM build from entry", true);
		zmap m1(entries);
		// t2.total();
		// auto vals = zmap::values(m1);
		// return std::make_tuple(m1, SFC_time);
		return m1;
	}	

	template<typename PT, typename M>
	auto map_insert(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();

		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();

		});
		
		parlay::sequence<par> insert_pts(n);
		parlay::parallel_for(0, n, [&](int i){
			insert_pts[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		auto m2 = zmap::multi_insert(mmp, insert_pts);

		return m2;
	}

	template<typename PT, typename M>
	auto map_delete(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();

		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});
		parlay::sequence<par> delete_pts(n);
		// parlay::sequence<pair<unsigned long long, long long> > delete_pts(n);
		parlay::parallel_for(0, n, [&](int i){
			delete_pts[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// delete_pts[i] = {P[i]->morton_id, P[i]->id};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		auto m2 = zmap::multi_delete(mmp, delete_pts);

		return m2;
	}

	template<class T, class MBR>
	auto range_count(T &zCPAM, MBR &query_mbr, bool use_hilbert = false){
		auto f = [&](auto &cur){ 
			return mbr_mbr_relation(cur, query_mbr);
		};

		auto f2 = [&](auto &cur){ 
			return point_in_mbr(cur, query_mbr);
		};

		// auto res = zmap::range_count_filter(zCPAM, f, f2);
		auto res = zmap::range_count_filter2(zCPAM, f, f2);
		return res;
		// auto ret = 0;
		// auto f2 = [&](auto cur){ 
		// 	auto flag = mbr_mbr_relation(cur.first, query_mbr);
		// 	if (flag < 0) {
		// 		return false;
		// 	}
		// 	else if (flag > 0) {
		// 		ret += cur.second;
		// 		return false;
		// 	}
		// 	else{
		// 		return true;
		// 	}
		// };
		// auto r2_map = zmap::aug_filter_mid(zCPAM, f2);

		// auto BL = query_mbr.first;
		// auto UR = query_mbr.second;
		// auto f3 = [&](auto cur){
		// 	if (BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
		// 		BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y){
		// 		return true;
		// 	}
		// 	return false;
		// };
		// r2_map = zmap::filter(r2_map, f3);
		// return ret + r2_map.size();
	}

	template<class T, class MBR>
	auto range_report(T &tree, MBR query_mbr, parlay::sequence<Point> &out, bool use_hilbert = false){
		// auto ret = zmap::values(filter_range(tree, query_mbr, use_hilbert));
		auto f = [&](auto cur){ 
			return mbr_mbr_relation(cur, query_mbr);
		};

		int64_t ret = 0;
		// zmap::range_report_filter(tree, f, ret, out);
		zmap::range_report_filter2(tree, f, ret, out);
		return ret;
	}

	//	return size of interior nodes and sizeof leaf nodes size, respectively
	auto size_in_bytes(){
		size_t inte_used = zmap::GC::used_node();
		size_t internal_nodes_space = sizeof(typename zmap::GC::regular_node) * inte_used;
		auto [used, unused] = parlay::internal::get_default_allocator().stats();
		return make_tuple(internal_nodes_space, used);
	}
}
