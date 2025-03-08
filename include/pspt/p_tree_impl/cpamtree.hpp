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

	template<typename M>
	auto map_diff(M &lhs, M &rhs){
		auto add = zmap::values(zmap::map_difference(rhs, lhs));
		auto remove = zmap::values(zmap::map_difference(lhs, rhs));
		return make_tuple(add, remove);
	}

	template<class PT>
	auto map_init(PT &P, bool use_hilbert = false){
		size_t n = P.size();
		
		parlay::parallel_for(0, n, [&](int i){
			P[i].morton_id = use_hilbert ? P[i].overlap_bits() : P[i].interleave_bits();
		});

		parlay::sequence<par> entries(n);
		parlay::parallel_for(0, n, [&](int i){
			entries[i] = {{P[i].morton_id, P[i].id}, P[i]};
			// entries[i] = {P[i]->id, P[i]};
		});
		zmap m1(entries);
		auto vals = zmap::values(m1);
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
	auto range_count(T &zCPAM, MBR query_mbr, bool use_hilbert = false){
		auto f = [&](auto cur, auto debug){ 
			if (debug == 1){
				cout << "cur mbr: "; print_mbr(cur);
				cout << "query mbr: "; print_mbr(query_mbr);
			}
			return mbr_mbr_relation(cur, query_mbr);
		};

		auto res = zmap::range_count_filter(zCPAM, f);
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
	auto range_report(T &zCPAM, MBR query_mbr, bool use_hilbert = false){
		auto BL = query_mbr.first;
		auto UR = query_mbr.second;
		// unsigned long long Z_min = 0, Z_max = 0;
		// if (!use_hilbert){
		// 	auto BR = geobase::Point(UR.x, BL.y), UL = geobase::Point(BL.x, UR.y);
		// 	Z_min = std::min(BL.interleave_bits(), std::min(std::min(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
		// 	Z_max = std::max(BL.interleave_bits(), std::max(std::max(UR.interleave_bits(), UL.interleave_bits()), BR.interleave_bits()));
		// }
		// else{
		// 	unsigned long long p1[] = {static_cast<unsigned long long>(BL.x), static_cast<unsigned long long>(BL.y)};	//	BL
		// 	unsigned long long p2[] = {static_cast<unsigned long long>(UR.x), static_cast<unsigned long long>(UR.y)};	//	UR
		// 	unsigned long long pointlo[] = {p1[0], p1[1]};
		// 	unsigned long long work[] = {p2[0], p2[1]};
  		// 	hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 1, pointlo, work);
		// 	Z_min = hilbert_c2i(2, 32, pointlo);
		// 	unsigned long long pointhi[] = {p2[0], p2[1]};
		// 	work[0] = p1[0], work[1] = p1[1];
  		// 	hilbert_box_pt(2, sizeof(unsigned long long), 8 * sizeof(unsigned long long), 0, work, pointhi);
		// 	Z_max = hilbert_c2i(2, 32, pointhi);
		// }
		// entry::key_t small(Z_min, 0);
		// entry::key_t large(Z_max, numeric_limits<unsigned long long>::max());
		// T r2_map = zmap::range(zCPAM, Z_min, Z_max);
		// T r2_map = zmap::range(zCPAM, small, large);
		// cout << r2_map.size() << endl;
		auto f2 = [&](auto cur){ 
			// print_mbr(cur);
			return !mbr_exclude_mbr(query_mbr, cur.first);
			// return !(UR.x < get<1>(cur).x || get<1>(cur).x < BL.x ||
			// 	UR.y < get<1>(cur).y || get<1>(cur).y < BL.y);
		};
		auto r2_map = zmap::aug_filter(zCPAM, f2);
		// r2_map = zmap::aug_filter(zCPAM, f2);
		// r2_map = zmap::aug_filter(r2_map, f2);

		auto f3 = [&](par cur){ 
			if (BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
				BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y){
				// cout << get<1>(cur).x << ", " << get<1>(cur).y << endl;
				return true;
			}
			return false;
		};

		r2_map = zmap::filter(r2_map, f3);

		auto ret = zmap::values(r2_map); 
		return ret;
	}

}
