#pragma once

#include "geobase.h"
#include "pspt/base_tree.h"
// #define USE_MBR
// #define SEQ
#define USE_PT

namespace ZD {

extern geobase::Bounding_Box largest_mbr;
extern geobase::break_down zd_build_break_down;
extern size_t maxSize;
extern double zd_leaf_copy_time;
extern double zd_inte_copy_time;

using namespace std;
using namespace geobase;
using parlay::par_do;
using parlay::par_do_if;
using parlay::sequence;

struct BaseNode {
#ifdef USE_MBR
  Bounding_Box mbr;
  BaseNode() : mbr({Point(-1, -1), Point(-1, -1)}) {}
#endif
  BaseNode() {}

  virtual ~BaseNode() = default;
  virtual bool is_leaf() { return false; }
  virtual size_t get_num_points() { return 0; }
};

struct InteNode : BaseNode {
  shared_ptr<BaseNode> l_son, r_son;
  size_t num_pts;

  InteNode() : l_son(nullptr), r_son(nullptr), num_pts(0) {}
  // InteNode(InteNode &x): l_son(x->l_son), r_son(x->r_son),
  // num_pts(x->num_pts){}

  virtual bool is_leaf() { return false; }
  virtual size_t get_num_points() { return num_pts; }
};

struct LeafNode : BaseNode {
  // sequence<Point> records;
  sequence<Point> records = sequence<Point>::uninitialized(32);

  template <typename Records>
  LeafNode(Records& r) {
    if (r.size() > 32) {
      records = sequence<Point>::uninitialized(r.size());
    }
    size_t i = 0;
    for (auto& pt : r) {
      parlay::assign_uninitialized(records[i++], pt);
      // records[i] = r[i];
    }
    records.resize(r.size());
  }

  template <typename Records, typename Func>
  LeafNode(Records& r, Func& f) {
    if (r.size() > 32) {
      records = sequence<Point>::uninitialized(r.size());
    }
    size_t i = 0;
    for (auto& pt : r) {
      if (f) {
        parlay::assign_uninitialized(records[i++], pt);
      }
      // records[i] = r[i];
    }
    records.resize(i);
  }

  // LeafNode(LeafNode &x): records(x->records){}

  virtual bool is_leaf() { return true; }
  virtual size_t get_num_points() { return records.size(); }

  void print_records() {
    cout << records.size() << endl;
    for (size_t i = 0; i < records.size(); i++) {
      cout << "(" << records[i].x << ", " << records[i].y << ")" << endl;
    }
  }
};

class Tree {
 public:
  size_t granularity_cutoff = 1000;
  size_t leaf_size = 1;

  shared_ptr<BaseNode> root;
  vector<shared_ptr<BaseNode> > multi_version_roots = {};

  Tree(size_t _leaf_sz);

  // entrance of building zdtree & tree construction
  shared_ptr<BaseNode> build(sequence<Point>& P, size_t l, size_t r, size_t b);
  void build(sequence<Point>& P);

  void clear();

  void merge_nodes(
      shared_ptr<BaseNode>& lhs, shared_ptr<BaseNode>& rhs,
      shared_ptr<InteNode>& cur);  // merge two sons to current node.
  void delete_merge_nodes(shared_ptr<BaseNode>& L, shared_ptr<BaseNode>& R,
                          InteNode* cur_node);

  shared_ptr<InteNode> create_internal(
      shared_ptr<BaseNode>& L,
      shared_ptr<BaseNode>& R);  //	create an internal node, do not store
                                 // pointers to original records.
  shared_ptr<LeafNode> create_leaf(
      sequence<Point>& P, size_t l, size_t r,
      size_t b);  // create a leaf, store all (pointers of) records.

  // 	in-place insertion
  void batch_insert_sorted(sequence<Point>& P);
  void batch_insert_sorted_node(shared_ptr<BaseNode>& x, sequence<Point>& P,
                                size_t l, size_t r, size_t b);

  //	in-place deletion
  void batch_delete_sorted(sequence<Point>& P);
  void batch_delete_sorted_node(shared_ptr<BaseNode>& x, sequence<Point>& P,
                                size_t l, size_t r, size_t b);

  // range report
  template <class Out>
  void range_report_node(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                         Bounding_Box& cur_mbr, FT x_prefix, FT y_prefix,
                         size_t b, bool x_splitter, size_t& cnt, Out& out);
  template <class Out>
  void range_report(Bounding_Box& query_mbr, Bounding_Box& cur_mbr, size_t& cnt,
                    Out& out);
  template <class Out>
  void range_report(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                    Bounding_Box& cur_mbr, size_t& cnt, Out& out);

  // range count
  size_t range_count_node(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                          Bounding_Box& cur_mbr, FT x_prefix, FT y_prefix,
                          size_t b, bool x_splitter);
  size_t range_count(Bounding_Box& query_mbr, Bounding_Box& cur_mbr);
  size_t range_count(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                     Bounding_Box& cur_mbr);

  // k nearest neighbor report
  template <class T>
  void knn_report_node(shared_ptr<BaseNode>& x, size_t& k, Point query_point,
                       Bounding_Box& cur_mbr, FT x_prefix, FT y_prefix,
                       size_t b, bool x_splitter, T& nn_res);
  auto knn_report(size_t& k, Point query_point, Bounding_Box& cur_mbr);
};

Tree::Tree(size_t _leaf_sz) { leaf_size = _leaf_sz; }

void Tree::clear() {
  root.reset();
  multi_version_roots.clear();
}

void Tree::delete_merge_nodes(shared_ptr<BaseNode>& L, shared_ptr<BaseNode>& R,
                              InteNode* cur_node) {
  // deal with MBR, covered points of parent
  auto L_num_pts = (L == nullptr) ? 0 : L->get_num_points();
  auto R_num_pts = (R == nullptr) ? 0 : R->get_num_points();

  cur_node->num_pts = L_num_pts + R_num_pts;
  cur_node->l_son = move(L);
  cur_node->r_son = move(R);
}

void Tree::merge_nodes(shared_ptr<BaseNode>& L, shared_ptr<BaseNode>& R,
                       shared_ptr<InteNode>& cur_node) {
  // deal with MBR, covered points of parent
  auto L_num_pts = L == nullptr ? 0 : L->get_num_points();
  auto R_num_pts = R == nullptr ? 0 : R->get_num_points();

  cur_node->num_pts = L_num_pts + R_num_pts;
  cur_node->l_son = move(L);
  cur_node->r_son = move(R);
}

shared_ptr<InteNode> Tree::create_internal(shared_ptr<BaseNode>& L,
                                           shared_ptr<BaseNode>& R) {
  shared_ptr<InteNode> cur_node(new InteNode());
  // augmented changes happen here
  merge_nodes(L, R, cur_node);
  return cur_node;
}

shared_ptr<LeafNode> Tree::create_leaf(sequence<Point>& P, size_t l, size_t r,
                                       size_t b) {
  auto cur_records = parlay::make_slice(&P[l], &P[r]);
  shared_ptr<LeafNode> cur_node(new LeafNode(cur_records));
  return cur_node;
}

shared_ptr<BaseNode> Tree::build(sequence<Point>& P, size_t l, size_t r,
                                 size_t b) {
  if (!b || (r - l <= leaf_size)) {
    return create_leaf(P, l, r, b);
  }

  auto splitter = split_by_bit(P, l, r, b);
  shared_ptr<BaseNode> L = nullptr;
  shared_ptr<BaseNode> R = nullptr;
  auto build_left = [&]() {
    if (l < splitter) L = build(P, l, splitter, b - 1);
  };
  auto build_right = [&]() {
    if (splitter < r) R = build(P, splitter, r, b - 1);
  };
  par_do_if(r - l >= granularity_cutoff, build_left, build_right);
  return create_internal(L, R);
}

void Tree::build(sequence<Point>& P) {
  if (!P.size()) return;
  root = build(P, 0, P.size(), 64);
}

void Tree::batch_insert_sorted_node(shared_ptr<BaseNode>& x, sequence<Point>& P,
                                    size_t l, size_t r, size_t b) {
  if (x == nullptr) {
    x = build(P, l, r, b);
    return;
  }
  auto less = [&](auto lhs, auto rhs) {
    return lhs.morton_id < rhs.morton_id ||
           (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
  };
  if (x->is_leaf()) {
    auto cur_leaf = static_cast<LeafNode*>(x.get());
    auto cur_records = parlay::make_slice(&P[l], &P[r]);
    if (!b || cur_leaf->records.size() + cur_records.size() <=
                  leaf_size) {  // current leaf is not full
      cur_leaf->records = parlay::merge(cur_leaf->records,
                                        parlay::make_slice(&P[l], &P[r]), less);
      return;
    } else {
      auto new_points = parlay::merge(cur_leaf->records,
                                      parlay::make_slice(&P[l], &P[r]), less);
      x = build(new_points, 0, new_points.size(), b);
      return;
    }
  }
  auto splitter = split_by_bit(P, l, r, b);
  auto cur_inte = static_cast<InteNode*>(x.get());
  auto insert_left = [&]() {
    if (l < splitter) {
      batch_insert_sorted_node(cur_inte->l_son, P, l, splitter, b - 1);
    };
  };
  auto insert_right = [&]() {
    if (splitter < r) {
      batch_insert_sorted_node(cur_inte->r_son, P, splitter, r, b - 1);
    };
  };
  par_do_if(r - l >= 256, insert_left, insert_right);
  delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
}

void Tree::batch_insert_sorted(sequence<Point>& P) {
  if (!P.size()) return;
  if (root == nullptr)
    build(P);
  else
    batch_insert_sorted_node(root, P, 0, P.size(), 64);
}

void Tree::batch_delete_sorted_node(shared_ptr<BaseNode>& x, sequence<Point>& P,
                                    size_t l, size_t r, size_t b) {
  if (x == nullptr) {
    return;
  }
  if (x->is_leaf()) {
    auto cur_leaf = static_cast<LeafNode*>(x.get());
    cur_leaf->records = get_delete_p(cur_leaf->records, P, l, r);
    if (!cur_leaf->records.size()) x.reset();
    return;
  }

  auto splitter = split_by_bit(P, l, r, b);
  auto cur_inte = static_cast<InteNode*>(x.get());
  auto delete_left = [&]() {
    if (l < splitter) {
      batch_delete_sorted_node(cur_inte->l_son, P, l, splitter, b - 1);
    };
  };
  auto delete_right = [&]() {
    if (splitter < r) {
      batch_delete_sorted_node(cur_inte->r_son, P, splitter, r, b - 1);
    };
  };

  par_do_if(r - l >= 256, delete_left, delete_right);

  auto less = [&](auto lhs, auto rhs) {
    return lhs.morton_id < rhs.morton_id ||
           (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
  };

  if (!cur_inte->l_son && !cur_inte->r_son)
    x.reset();
  else {
    if (!cur_inte->l_son) {
      if (cur_inte->r_son->get_num_points() <= leaf_size)
        x = move(cur_inte->r_son);
    } else {
      if (!cur_inte->r_son) {
        if (cur_inte->l_son->get_num_points() <= leaf_size)
          x = move(cur_inte->l_son);
      } else {
        if (cur_inte->l_son->get_num_points() +
                cur_inte->r_son->get_num_points() <=
            leaf_size) {
          auto L = static_cast<LeafNode*>(cur_inte->l_son.get());
          auto R = static_cast<LeafNode*>(cur_inte->r_son.get());
          auto cur_records = parlay::merge(L->records, R->records, less);
          x = create_leaf(cur_records, 0, cur_records.size(), 0);
        } else {
          delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
        }
      }
    }
    auto L_num_pts =
        cur_inte->l_son == nullptr ? 0 : cur_inte->l_son->get_num_points();
    auto R_num_pts =
        cur_inte->r_son == nullptr ? 0 : cur_inte->r_son->get_num_points();
    cur_inte->num_pts = L_num_pts + R_num_pts;
  }
}

void Tree::batch_delete_sorted(sequence<Point>& P) {
  if (!P.size() || root == nullptr)
    return;
  else
    batch_delete_sorted_node(root, P, 0, P.size(), 64);
}

size_t Tree::range_count_node(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                              Bounding_Box& cur_mbr, FT x_prefix, FT y_prefix,
                              size_t b, bool x_splitter) {
  int flag = mbr_mbr_relation(cur_mbr, query_mbr);
  if (flag < 0) return 0;
  if (flag > 0) {
    return x->get_num_points();
  }
  if (x->is_leaf()) {  // we have to scan the leaf to report the number of
                       // points;
    size_t ret = 0;
    auto cur_leaf = static_cast<LeafNode*>(x.get());
    for (auto& p : cur_leaf->records) {
      if (point_in_mbr(p, query_mbr)) {
        ret += 1;
      }
    }
    return ret;
  } else {
    auto cur_inte = static_cast<InteNode*>(x.get());
    FT splitter = 1.0 * (1u << (b - 1));
    size_t ret_L = 0, ret_R = 0;
    if (cur_inte->l_son != nullptr) {
      auto L_box = cur_mbr;
      if (x_splitter) {
        L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
        ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix,
                                 y_prefix, b, !x_splitter);
      } else {
        L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
        ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix,
                                 y_prefix, b - 1, !x_splitter);
      }
    }
    if (cur_inte->r_son != nullptr) {
      auto R_box = cur_mbr;
      if (x_splitter) {
        R_box.first.x = max(x_prefix + splitter, R_box.first.x);
        ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box,
                                 x_prefix + splitter, y_prefix, b, !x_splitter);
      } else {
        R_box.first.y = max(y_prefix + splitter, R_box.first.y);
        ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, x_prefix,
                                 y_prefix + splitter, b - 1, !x_splitter);
      }
    }
    return ret_L + ret_R;
  }
  return -1;  // unexpected error happens if the code runs to here.
}

size_t Tree::range_count(Bounding_Box& query_mbr, Bounding_Box& cur_mbr) {
  size_t ret = range_count_node(root, query_mbr, cur_mbr, 0.0, 0.0, 32, true);
  return ret;
}

size_t Tree::range_count(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                         Bounding_Box& cur_mbr) {
  size_t ret = range_count_node(x, query_mbr, cur_mbr, 0.0, 0.0, 32, true);
  return ret;
}

template <class Out>
void Tree::range_report_node(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                             Bounding_Box& cur_mbr, FT x_prefix, FT y_prefix,
                             size_t b, bool x_splitter, size_t& cnt, Out& out) {
  if (!x) {
    return;
  }
  auto flag = mbr_mbr_relation(cur_mbr, query_mbr);
  if (flag < 0) return;

  if (x->is_leaf()) {
    auto cur_leaf = static_cast<LeafNode*>(x.get());
    for (auto& p : cur_leaf->records) {
      if (point_in_mbr(p, query_mbr)) {
        out[cnt++] = p;
      }
    }
    return;
  }
  auto cur_inte = static_cast<InteNode*>(x.get());
  auto [L_box, R_box, rx_prefix, ry_prefix] =
      compute_cur_box(cur_mbr, x_prefix, y_prefix, b, x_splitter);
  range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix,
                    b - 1, !x_splitter, cnt, out);
  range_report_node(cur_inte->r_son, query_mbr, R_box, rx_prefix, ry_prefix,
                    b - 1, !x_splitter, cnt, out);
}

template <class Out>
void Tree::range_report(Bounding_Box& query_mbr, Bounding_Box& cur_mbr,
                        size_t& cnt, Out& out) {
  range_report_node(root, query_mbr, cur_mbr, 0.0, 0.0, 64, true, cnt, out);
}

template <class Out>
void Tree::range_report(shared_ptr<BaseNode>& x, Bounding_Box& query_mbr,
                        Bounding_Box& cur_mbr, size_t& cnt, Out& out) {
  range_report_node(x, query_mbr, cur_mbr, 0.0, 0.0, 32, true, cnt, out);
}

auto Tree::knn_report(size_t& k, Point query_point, Bounding_Box& cur_mbr) {
  priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
  knn_report_node(root, k, query_point, cur_mbr, 0.0, 0.0, 32, true, nn_res);
  return nn_res;
}

template <class T>
void Tree::knn_report_node(shared_ptr<BaseNode>& x, size_t& k,
                           Point query_point, Bounding_Box& cur_mbr,
                           FT x_prefix, FT y_prefix, size_t b, bool x_splitter,
                           T& nn_res) {
  if (!x) return;
  if (x->is_leaf()) {
    auto cur_leaf = static_cast<LeafNode*>(x.get());
    for (auto& p : cur_leaf->records) {
      auto cur_sqrdis = point_point_sqrdis(p, query_point);
      if (nn_res.size() < k) {
        nn_res.push({p, cur_sqrdis});
      } else if (cur_sqrdis < nn_res.top().second) {
        nn_res.pop();
        nn_res.push({p, cur_sqrdis});
      }
    }
    return;
  }
  auto cur_inte = static_cast<InteNode*>(x.get());
  auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
  FT splitter = 1.0 * (1u << (b - 1));
  auto L_box = cur_mbr, R_box = cur_mbr;
  if (cur_inte->l_son != nullptr) {
    if (x_splitter) {
      L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
    } else {
      L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
    }
    l_son_sqrdis = point_mbr_sqrdis(query_point, L_box);
  }
  if (cur_inte->r_son != nullptr) {
    if (x_splitter) {
      R_box.first.x = max(x_prefix + splitter, R_box.first.x);
    } else {
      R_box.first.y = max(y_prefix + splitter, R_box.first.y);
    }
    r_son_sqrdis = point_mbr_sqrdis(query_point, R_box);
  }

  if (l_son_sqrdis <= r_son_sqrdis) {  // first go left
    if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second) {
      if (x_splitter)
        knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix,
                        y_prefix, b, !x_splitter, nn_res);
      else
        knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix,
                        y_prefix, b - 1, !x_splitter, nn_res);
    }
    if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second) {
      if (x_splitter)
        knn_report_node(cur_inte->r_son, k, query_point, R_box,
                        x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
      else
        knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix,
                        y_prefix + splitter, b - 1, !x_splitter, nn_res);
    }
  } else {  // first go right
    if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second) {
      if (x_splitter)
        knn_report_node(cur_inte->r_son, k, query_point, R_box,
                        x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
      else
        knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix,
                        y_prefix + splitter, b - 1, !x_splitter, nn_res);
    }
    if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second) {
      if (x_splitter)
        knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix,
                        y_prefix, b, !x_splitter, nn_res);
      else
        knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix,
                        y_prefix, b - 1, !x_splitter, nn_res);
    }
  }
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight = 6,
          uint_fast8_t kImbaRatio = 30>
class Zdtree
    : public pspt::BaseTree<Point,
                            Zdtree<Point, SplitRule, kSkHeight, kImbaRatio>,
                            kSkHeight, kImbaRatio> {
 public:
  using BT =
      pspt::BaseTree<Point, Zdtree<Point, SplitRule, kSkHeight, kImbaRatio>,
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
  using Leaf = pspt::Node;
  using Interior = pspt::Node;

  void DeleteTree() {
    tree.clear();
    return;
  }

  template <typename Range>
  void Build(Range In) {
    Slice A = parlay::make_slice(In);
    // cpam_aug_map_ = std::move(map_init(A));
  }
  void BatchInsert(Slice In) {
    // if (!this->cpam_aug_map_.root || !zmap::size(this->cpam_aug_map_.root)) {
    //   return Build(std::forward<Slice>(In));
    // }
    // cpam_aug_map_ = map_insert(In, std::move(cpam_aug_map_));
  }

  void BatchDelete(Slice In) {
    // if (!this->cpam_aug_map_.root || !zmap::size(this->cpam_aug_map_.root)) {
    //   return;
    // }
    // cpam_aug_map_ = map_delete(In, std::move(cpam_aug_map_));
  }

  template <typename Node, typename Range>
  auto KNN(Node* T, Point const& q, pspt::kBoundedQueue<Point, Range>& bq) {
    pspt::KNNLogger logger;
    // this->knn(this->cpam_aug_map_, q, bq.max_size(), logger.vis_leaf_num);
    // zmap::template knn<BT>(this->cpam_aug_map_, q, bq, logger);
    return logger;
  }

  auto RangeCount(Box const& q) {
    pspt::RangeQueryLogger logger;
    int size = 0;
    // auto size = this->range_count(this->cpam_aug_map_, q);
    return std::make_pair(size, logger);
  }

  template <typename Range>
  auto RangeQuery(Box const& q, Range&& Out) {
    pspt::RangeQueryLogger logger;
    // auto size = this->range_report(this->cpam_aug_map_, q, Out);
    int size = 0;
    return std::make_pair(size, logger);
  }

  constexpr static char const* GetTreeName() { return "ZdTree"; }
  constexpr static char const* CheckHasBox() { return "WithoutBox"; }
  Tree tree = Tree(32);
};

}  // namespace ZD
