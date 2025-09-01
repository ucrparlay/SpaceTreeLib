#pragma once

#include <cstdint>

#include "../psi/base_tree.h"
#include "cpam/cpam.h"
#include "dependence/loggers.h"
// #include "geobase.h"
// #include "hilbert.h"
#include "psi/base_tree.h"

namespace CPAMTree {
using namespace std;
// using namespace geobase;
using parlay::par_do;
using parlay::par_do_if;
using parlay::sequence;

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class CpamRaw
    : public psi::BaseTree<Point,
                            CpamRaw<Point, SplitRule, kSkHeight, kImbaRatio>,
                            kSkHeight, kImbaRatio> {
 public:
  using BT =
      psi::BaseTree<Point, CpamRaw<Point, SplitRule, kSkHeight, kImbaRatio>,
                     kSkHeight, kImbaRatio>;

  using BucketType = typename BT::BucketType;
  using BallsType = typename BT::BallsType;
  using DimsType = typename BT::DimsType;
  using BucketSeq = typename BT::BucketSeq;
  using BallSeq = typename BT::BallSeq;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using Slice = typename BT::Slice;
  using Points = typename BT::Points;
  using PointsIter = typename BT::PointsIter;
  using Box = typename BT::Box;
  using BoxSeq = typename BT::BoxSeq;
  using Circle = typename BT::NormalCircle;

  using HyperPlane = typename BT::HyperPlane;
  using HyperPlaneSeq = typename BT::HyperPlaneSeq;
  using NodeTag = typename BT::NodeTag;
  using NodeTagSeq = typename BT::NodeTagSeq;
  using NodeBox = typename BT::NodeBox;
  using NodeBoxSeq = typename BT::NodeBoxSeq;

  // begin of CPAM_RAW definition
  using key_type =
      pair<unsigned long long, unsigned long long>;  // morton_id, id
  using val_type = Point;
  using aug_type = pair<Box, size_t>;

  typedef pair<Point, int64_t> nn_pair;

  struct nn_pair_cmp {
    bool operator()(nn_pair& lhs, nn_pair& rhs) {
      return lhs.second < rhs.second ||
             (lhs.second == rhs.second && lhs.first.aug.id > rhs.first.aug.id);
    }
  };

  //	CPAM entry
  struct entry {
    using key_t = key_type;
    using val_t = val_type;
    using aug_t = aug_type;

    static inline bool comp(key_t a, key_t b) { return a < b; }
    static aug_t get_empty() { return make_pair(BT::GetEmptyBox(), 0); }
    static aug_t from_entry(key_t k, val_t v) {
      return make_pair(Box(v, v), 1);
    }
    static aug_t combine(aug_t a, aug_t b) {
      return make_pair(BT::GetBox(a.first, b.first), a.second + b.second);
    }
  };

  using zmap = cpam::aug_map<entry, 40>;
  using par = std::tuple<typename entry::key_t, typename entry::val_t>;

  using Leaf = zmap::Tree::node;
  using Interior = zmap::Tree::node;

  // template <typename T>
  // auto knn(T& tree, Point const& query_point, auto k, auto& vis_leaf) {
  //   auto f = [&](auto cur_pt) {
  //     return BT::P2PDistanceSquare(cur_pt, query_point);
  //     // return point_point_sqrdis(cur_pt, query_point);
  //   };
  //
  //   auto f2 = [&](auto cur_mbr) {
  //     return BT::P2BMinDistanceSquare(query_point, cur_mbr);
  //     // return point_mbr_sqrdis(query_point, cur_mbr);
  //   };
  //
  //   std::priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
  //   zmap::knn_filter(tree, f, f2, k, nn_res, vis_leaf);
  //   return nn_res.top().second;
  // }

  template <typename M>
  auto map_diff(M& lhs, M& rhs) {
    auto add = zmap::values(zmap::map_difference(rhs, lhs));
    auto remove = zmap::values(zmap::map_difference(lhs, rhs));
    return make_tuple(add, remove);
  }

  template <class PT>
  auto map_init(PT& P) {
    size_t n = P.size();

    // parlay::internal::timer t("SFC time", true);
    parlay::parallel_for(
        0, n, [&](int i) { P[i].aug.code = SplitRule::Encode(P[i]); });
    // t.stop();
    // t.total();
    // auto SFC_time = t.total_time();
    // auto h_values = parlay::sequence<unsigned long long>::uninitialized(n);
    // auto h_values2 = parlay::sequence<unsigned long long>::uninitialized(n);
    // parlay::internal::timer t1("Create entry time", true);
    parlay::sequence<par> entries(n);
    parlay::parallel_for(0, n, [&](int i) {
      entries[i] = {{P[i].aug.code, P[i].aug.id}, P[i]};
      // entries[i] = {P[i]->id, P[i]};
      // h_values[i] = P[i].morton_id;
      // h_values2[i] = P[i].morton_id;
    });
    // t1.total();

    // parlay::internal::timer t_integer_sort("integer_sort", true);
    // auto ret = parlay::integer_sort(h_values);
    // t_integer_sort.total();
    // auto less = [&](int a, int b) { return a < b; };
    // parlay::internal::timer t_sample_sort("sample_sort", true);
    // auto ret2 = parlay::internal::sample_sort(
    //     parlay::make_slice(h_values2.begin(), h_values2.end()), less);
    // t_sample_sort.total();

    // parlay::internal::timer t2("CPAM build from entry", true);
    zmap m1(entries);
    // t2.total();
    // auto vals = zmap::values(m1);
    // return std::make_tuple(m1, SFC_time);
    return std::move(m1);
  }

  template <typename PT, typename M>
  auto map_insert(PT P, M&& mmp) {
    size_t n = P.size();

    parlay::parallel_for(
        0, n, [&](int i) { P[i].aug.code = SplitRule::Encode(P[i]); });

    parlay::sequence<par> insert_pts(n);
    parlay::parallel_for(0, n, [&](int i) {
      insert_pts[i] = {{P[i].aug.code, P[i].aug.id}, P[i]};
      // insert_pts[i] = {P[i]->id, P[i]};
    });

    return std::move(zmap::multi_insert(std::forward<M>(mmp), insert_pts));
  }

  template <typename PT, typename M>
  auto map_delete(PT P, M&& mmp) {
    size_t n = P.size();

    parlay::parallel_for(
        0, n, [&](int i) { P[i].aug.code = SplitRule::Encode(P[i]); });
    parlay::sequence<par> delete_pts(n);
    // parlay::sequence<pair<unsigned long long, long long> > delete_pts(n);
    parlay::parallel_for(0, n, [&](int i) {
      delete_pts[i] = {{P[i].aug.code, P[i].aug.id}, P[i]};
      // delete_pts[i] = {P[i]->morton_id, P[i]->id};
      // insert_pts[i] = {P[i]->id, P[i]};
    });

    return std::move(zmap::multi_delete(std::forward<M>(mmp), delete_pts));
  }

  int mbr_mbr_relation(Box const& a, Box const& b) {
    if (!BT::BoxIntersectBox(a, b)) {
      return -1;  // disjoint
    } else if (BT::WithinBox(a, b)) {
      return 1;
    } else {
      return 0;  // overlap
    }
  }

  template <class T, class MBR>
  auto range_count(T& zCPAM, MBR& query_mbr) {
    // auto f = [&](auto& cur) { return mbr_mbr_relation(cur, query_mbr); };
    // auto f2 = [&](auto& cur) { return point_in_mbr(cur, query_mbr); };
    auto f = [&](auto& cur_mbr) {
      return mbr_mbr_relation(cur_mbr, query_mbr);
    };
    auto f2 = [&](auto& cur_point) {
      return BT::WithinBox(cur_point, query_mbr);
    };

    auto res = zmap::range_count_filter2(zCPAM, f, f2);
    return res;
  }

  template <class T, class MBR, typename Out>
  auto range_report(T& tree, MBR query_mbr, Out&& out) {
    // auto ret = zmap::values(filter_range(tree, query_mbr, use_hilbert));
    auto f = [&](auto cur) { return mbr_mbr_relation(cur, query_mbr); };

    int64_t ret = 0;
    // zmap::range_report_filter(tree, f, ret, out);
    zmap::range_report_filter2(tree, f, ret, out);
    return ret;
  }

  //	return size of interior nodes and sizeof leaf nodes size, respectively
  auto size_in_bytes() {
    size_t inte_used = zmap::GC::used_node();
    size_t internal_nodes_space =
        sizeof(typename zmap::GC::regular_node) * inte_used;
    auto [used, unused] = parlay::internal::get_default_allocator().stats();
    return make_tuple(internal_nodes_space, used);
  }

  // APIs
  template <typename Range>
  void Build(Range In) {
    Slice A = parlay::make_slice(In);
    cpam_aug_map_ = std::move(map_init(A));
  }

  void BatchInsert(Slice In) {
    if (!this->cpam_aug_map_.root || !zmap::size(this->cpam_aug_map_.root)) {
      return Build(std::forward<Slice>(In));
    }
    cpam_aug_map_ = map_insert(In, std::move(cpam_aug_map_));
  }

  void BatchDelete(Slice In) {
    if (!this->cpam_aug_map_.root || !zmap::size(this->cpam_aug_map_.root)) {
      return;
    }
    cpam_aug_map_ = map_delete(In, std::move(cpam_aug_map_));
  }

  template <typename Node, typename Range>
  auto KNN(Node* T, Point const& q, psi::kBoundedQueue<Point, Range>& bq) {
    psi::KNNLogger logger;
    // this->knn(this->cpam_aug_map_, q, bq.max_size(), logger.vis_leaf_num);
    zmap::template knn<BT>(this->cpam_aug_map_, q, bq, logger);
    return logger;
  }

  auto RangeCount(Box const& q) {
    psi::RangeQueryLogger logger;
    auto size = this->range_count(this->cpam_aug_map_, q);
    return std::make_pair(size, logger);
  }

  template <typename Range>
  auto RangeQuery(Box const& q, Range&& Out) {
    psi::RangeQueryLogger logger;
    auto size = this->range_report(this->cpam_aug_map_, q, Out);
    return std::make_pair(size, logger);
  }

  void DeleteTree() {
    cpam_aug_map_.clear();
    return;
  }

  size_t GetSize() const { return cpam_aug_map_.size(); }
  constexpr static char const* GetTreeName() { return "CPAM"; }
  constexpr static char const* CheckHasBox() { return "HasBox"; }
  zmap cpam_aug_map_;
};

}  // namespace CPAMTree
