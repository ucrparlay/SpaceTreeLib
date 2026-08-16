#pragma once

#include "geobase_3d.h"
#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"
#define USE_MBR
// #define SEQ
// #define USE_PT

namespace ZD3D
{

extern geobase::bounding_box_type largest_mbr;
extern size_t maxSize;
extern double zd_leaf_copy_time;
extern double zd_inte_copy_time;

using namespace std;
using namespace geobase;
using parlay::par_do;
using parlay::par_do_if;
using parlay::sequence;

bounding_box_type empty_mbr(Point(ft_inf_max, ft_inf_max, ft_inf_max),
			    Point(ft_inf_min, ft_inf_min, ft_inf_min));

struct base_node {
#ifdef USE_MBR
	bounding_box_type mbr;
	base_node() : mbr(empty_mbr)
	{
	}
#endif
	// base_node() {}

	virtual ~base_node() = default;
	virtual bool is_leaf()
	{
		return false;
	}
	virtual size_t get_num_points()
	{
		return 0;
	}
};

struct inte_node : base_node {
	shared_ptr<base_node> l_son, r_son;
	size_t num_pts;

	inte_node() : l_son(nullptr), r_son(nullptr), num_pts(0)
	{
	}
	// inte_node(inte_node &x): l_son(x->l_son), r_son(x->r_son),
	// num_pts(x->num_pts){}

	virtual bool is_leaf()
	{
		return false;
	}
	virtual size_t get_num_points()
	{
		return num_pts;
	}
};

struct leaf_node : base_node {
	// sequence<Point> records;
	sequence<Point> records = sequence<Point>::uninitialized(32);

	template <typename Records>
	leaf_node(Records &r)
	{
		if (r.size() > 32) {
			records = sequence<Point>::uninitialized(r.size());
		}
		size_t i = 0;
		for (auto &pt : r) {
			parlay::assign_uninitialized(records[i++], pt);
			// records[i] = r[i];
		}
		records.resize(r.size());
		mbr = get_mbr(r);
	}

	template <typename Records, typename Func>
	leaf_node(Records &r, Func &f)
	{
		if (r.size() > 32) {
			records = sequence<Point>::uninitialized(r.size());
		}
		size_t i = 0;
		for (auto &pt : r) {
			if (f) {
				parlay::assign_uninitialized(records[i++], pt);
			}
			// records[i] = r[i];
		}
		records.resize(i);
		mbr = get_mbr(r);
	}

	// leaf_node(leaf_node &x): records(x->records){}

	virtual bool is_leaf()
	{
		return true;
	}
	virtual size_t get_num_points()
	{
		return records.size();
	}

	void print_records()
	{
		cout << records.size() << endl;
		for (size_t i = 0; i < records.size(); i++) {
			cout << "(" << records[i].x << ", " << records[i].y
			     << ")" << endl;
		}
	}
};

class Tree
{
public:
	size_t granularity_cutoff = 1000;
	size_t leaf_size = 1;

	shared_ptr<base_node> root;
	vector<shared_ptr<base_node>> multi_version_roots = {};

	Tree(size_t _leaf_sz);

	// entrance of building zdtree & tree construction
	shared_ptr<base_node> build(sequence<Point> &P, size_t l, size_t r,
				    size_t b);
	void build(sequence<Point> &P);

	auto collect_records(shared_ptr<base_node> &x);

	void clear();

	void merge_nodes(
		shared_ptr<base_node> &lhs, shared_ptr<base_node> &rhs,
		shared_ptr<inte_node> &cur); // merge two sons to current node.
	void delete_merge_nodes(shared_ptr<base_node> &L,
				shared_ptr<base_node> &R, inte_node *cur_node);

	shared_ptr<inte_node>
	create_internal(shared_ptr<base_node> &L,
			shared_ptr<base_node>
				&R); //	create an internal node, do not store
				     // pointers to original records.
	shared_ptr<leaf_node> create_leaf(
		sequence<Point> &P, size_t l, size_t r,
		size_t b); // create a leaf, store all (pointers of) records.

	// 	in-place insertion
	void batch_insert_sorted(sequence<Point> &P);
	void batch_insert_sorted_node(shared_ptr<base_node> &x,
				      sequence<Point> &P, size_t l, size_t r,
				      size_t b);

	//	in-place deletion
	void batch_delete_sorted(sequence<Point> &P);
	void batch_delete_sorted_node(shared_ptr<base_node> &x,
				      sequence<Point> &P, size_t l, size_t r,
				      size_t b);

	// range report
	template <class Out>
	void range_report_node(shared_ptr<base_node> &x,
			       bounding_box_type &query_mbr, size_t &cnt,
			       Out &out);
	template <class Out>
	void range_report(bounding_box_type &query_mbr, size_t &cnt, Out &out);

	// range count
	size_t range_count_node(shared_ptr<base_node> &x,
				bounding_box_type &query_mbr);
	size_t range_count(bounding_box_type &query_mbr);

	// k nearest neighbor report
	template <class T>
	void knn_report_node(shared_ptr<base_node> &x, size_t &k,
			     Point query_point, T &nn_res);
	auto knn_report(size_t &k, Point query_point);

	auto check_mbr(shared_ptr<base_node> &x);
};

Tree::Tree(size_t _leaf_sz)
{
	leaf_size = _leaf_sz;
}

void Tree::clear()
{
	root.reset();
	multi_version_roots.clear();
}

void Tree::delete_merge_nodes(shared_ptr<base_node> &L,
			      shared_ptr<base_node> &R, inte_node *cur_node)
{
	// deal with MBR, covered points of parent
	auto l_num_pts = (L == nullptr) ? 0 : L->get_num_points();
	auto r_num_pts = (R == nullptr) ? 0 : R->get_num_points();
	auto l_mbr = L == nullptr ? empty_mbr : L->mbr;
	auto r_mbr = R == nullptr ? empty_mbr : R->mbr;

	cur_node->mbr = merge_mbr(l_mbr, r_mbr);
	cur_node->num_pts = l_num_pts + r_num_pts;
	cur_node->l_son = move(L);
	cur_node->r_son = move(R);
}

void Tree::merge_nodes(shared_ptr<base_node> &L, shared_ptr<base_node> &R,
		       shared_ptr<inte_node> &cur_node)
{
	// deal with MBR, covered points of parent
	auto l_num_pts = L == nullptr ? 0 : L->get_num_points();
	auto r_num_pts = R == nullptr ? 0 : R->get_num_points();
	auto l_mbr = L == nullptr ? empty_mbr : L->mbr;
	auto r_mbr = R == nullptr ? empty_mbr : R->mbr;

	cur_node->mbr = merge_mbr(l_mbr, r_mbr);
	cur_node->num_pts = l_num_pts + r_num_pts;
	cur_node->l_son = move(L);
	cur_node->r_son = move(R);
}

shared_ptr<inte_node> Tree::create_internal(shared_ptr<base_node> &L,
					    shared_ptr<base_node> &R)
{
	shared_ptr<inte_node> cur_node(new inte_node());
	// augmented changes happen here
	merge_nodes(L, R, cur_node);
	return cur_node;
}

shared_ptr<leaf_node> Tree::create_leaf(sequence<Point> &P, size_t l, size_t r,
					size_t b)
{
	auto cur_records = parlay::make_slice(&P[l], &P[r]);
	shared_ptr<leaf_node> cur_node(new leaf_node(cur_records));
	return cur_node;
}

shared_ptr<base_node> Tree::build(sequence<Point> &P, size_t l, size_t r,
				  size_t b)
{
	// cout << "l, r = " << l << ", " << r << endl;
	// cout << "b = " << b << endl;
	if (!b || (r - l <= leaf_size)) {
		// if (b == 1 || (r - l <= leaf_size)) {
		return create_leaf(P, l, r, b);
	}

	auto splitter = split_by_bit(P, l, r, b);
	// cout << "splitter = " << splitter << endl;
	shared_ptr<base_node> L = nullptr;
	shared_ptr<base_node> R = nullptr;
	auto build_left = [&]() {
		if (l < splitter)
			L = build(P, l, splitter, b - 1);
	};
	auto build_right = [&]() {
		if (splitter < r)
			R = build(P, splitter, r, b - 1);
	};
	par_do_if(r - l >= granularity_cutoff, build_left, build_right);
	return create_internal(L, R);
}

void Tree::build(sequence<Point> &P)
{
	if (!P.size())
		return;
	root = build(P, 0, P.size(), 63);
}

void Tree::batch_insert_sorted_node(shared_ptr<base_node> &x,
				    sequence<Point> &P, size_t l, size_t r,
				    size_t b)
{
	if (x == nullptr) {
		x = build(P, l, r, b);
		return;
	}
	auto less = [&](auto lhs, auto rhs) {
		return lhs.morton_id < rhs.morton_id ||
		       (lhs.morton_id == rhs.morton_id && lhs.id < rhs.id);
	};
	if (x->is_leaf()) {
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		auto cur_records = parlay::make_slice(&P[l], &P[r]);
		if (!b || cur_leaf->records.size() + cur_records.size() <=
				  leaf_size) { // current leaf is not full
			// if (b == 1 || cur_leaf->records.size() +
			// cur_records.size() <= leaf_size) {  // current leaf
			// is not full
			cur_leaf->records = parlay::merge(
				cur_leaf->records,
				parlay::make_slice(&P[l], &P[r]), less);
			cur_leaf->mbr = get_mbr(cur_leaf->records);
			return;
		} else {
			auto new_points = parlay::merge(
				cur_leaf->records,
				parlay::make_slice(&P[l], &P[r]), less);
			x = build(new_points, 0, new_points.size(), b);
			return;
		}
	}
	auto splitter = split_by_bit(P, l, r, b);
	auto cur_inte = static_cast<inte_node *>(x.get());
	auto insert_left = [&]() {
		if (l < splitter) {
			batch_insert_sorted_node(cur_inte->l_son, P, l,
						 splitter, b - 1);
		};
	};
	auto insert_right = [&]() {
		if (splitter < r) {
			batch_insert_sorted_node(cur_inte->r_son, P, splitter,
						 r, b - 1);
		};
	};
	par_do_if(r - l >= 256, insert_left, insert_right);
	delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
}

void Tree::batch_insert_sorted(sequence<Point> &P)
{
	if (!P.size())
		return;
	if (root == nullptr)
		build(P);
	else
		batch_insert_sorted_node(root, P, 0, P.size(), 63);
}

void Tree::batch_delete_sorted_node(shared_ptr<base_node> &x,
				    sequence<Point> &P, size_t l, size_t r,
				    size_t b)
{
	if (!x) {
		return;
	}
	if (x->is_leaf()) {
		// cout << "leaf_type-delete: " << r - l << " points, b = " << b
		// << endl;
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		// cout << "cur mbr: "; print_mbr(x->mbr);
		// cout << "before size: " << cur_leaf->get_num_points() <<
		// endl;
		cur_leaf->records = get_delete_p(cur_leaf->records, P, l, r);
		// cout << "after size: " << cur_leaf->get_num_points() << endl;
		cur_leaf->mbr = get_mbr(cur_leaf->records);
		if (!cur_leaf->records.size())
			x.reset();
		return;
	}

	// cout << "calculate splitter" << endl;
	auto splitter = split_by_bit(P, l, r, b);
	// cout << "splitter = " << splitter << endl;
	auto cur_inte = static_cast<inte_node *>(x.get());
	auto delete_left = [&]() {
		if (l < splitter) {
			batch_delete_sorted_node(cur_inte->l_son, P, l,
						 splitter, b - 1);
		};
	};
	auto delete_right = [&]() {
		if (splitter < r) {
			batch_delete_sorted_node(cur_inte->r_son, P, splitter,
						 r, b - 1);
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
				if (cur_inte->l_son->get_num_points() <=
				    leaf_size)
					x = move(cur_inte->l_son);
			} else {
				if (cur_inte->l_son->get_num_points() +
					    cur_inte->r_son->get_num_points() <=
				    leaf_size) {
					auto L = static_cast<leaf_node *>(
						cur_inte->l_son.get());
					auto R = static_cast<leaf_node *>(
						cur_inte->r_son.get());
					auto cur_records = parlay::merge(
						L->records, R->records, less);
					x = create_leaf(cur_records, 0,
							cur_records.size(), 0);
				} else {
					delete_merge_nodes(cur_inte->l_son,
							   cur_inte->r_son,
							   cur_inte);
				}
			}
		}
		auto l_num_pts = cur_inte->l_son == nullptr
					 ? 0
					 : cur_inte->l_son->get_num_points();
		auto l_mbr = cur_inte->l_son == nullptr ? empty_mbr
							: cur_inte->l_son->mbr;
		auto r_num_pts = cur_inte->r_son == nullptr
					 ? 0
					 : cur_inte->r_son->get_num_points();
		auto r_mbr = cur_inte->r_son == nullptr ? empty_mbr
							: cur_inte->r_son->mbr;

		cur_inte->mbr = merge_mbr(l_mbr, r_mbr);
		cur_inte->num_pts = l_num_pts + r_num_pts;
	}
}

void Tree::batch_delete_sorted(sequence<Point> &P)
{
	// cout << "delete " << P.size() << " points" << endl;
	if (!P.size() || root == nullptr)
		return;
	else
		batch_delete_sorted_node(root, P, 0, P.size(), 63);
}

size_t Tree::range_count_node(shared_ptr<base_node> &x,
			      bounding_box_type &query_mbr)
{
	// int flag = mbr_mbr_relation(x->mbr, query_mbr);
	int flag = zy_mbr_relation(x->mbr, query_mbr);
	if (flag < 0)
		return 0;
	if (flag > 0) {
		return x->get_num_points();
	}
	if (x->is_leaf()) { // we have to scan the leaf to report the number of
			    // points;
		size_t ret = 0;
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		for (auto &p : cur_leaf->records) {
			if (point_in_mbr(p, query_mbr)) {
				ret += 1;
			}
		}
		return ret;
	} else {
		auto cur_inte = static_cast<inte_node *>(x.get());
		size_t ret_L = 0, ret_R = 0;
		if (cur_inte->l_son != nullptr) {
			ret_L = range_count_node(cur_inte->l_son, query_mbr);
		}
		if (cur_inte->r_son != nullptr) {
			ret_R = range_count_node(cur_inte->r_son, query_mbr);
		}
		return ret_L + ret_R;
	}
	return -1; // unexpected error happens if the code runs to here.
}

size_t Tree::range_count(bounding_box_type &query_mbr)
{
	size_t ret = 0;
	ret = range_count_node(root, query_mbr);
	return ret;
}

template <class Out>
void Tree::range_report_node(shared_ptr<base_node> &x,
			     bounding_box_type &query_mbr, size_t &cnt,
			     Out &out)
{
	if (!x) {
		return;
	}
	// auto flag = mbr_mbr_relation(x->mbr, query_mbr);
	auto flag = zy_mbr_relation(x->mbr, query_mbr);
	if (flag < 0)
		return;

	if (x->is_leaf()) {
		// cout << "scanned leaf size: " << x->get_num_points() << endl;
		// if (x->get_num_points() > 32) {
		//   cout << "large leaf touched" << endl;
		// }
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		for (auto &p : cur_leaf->records) {
			if (point_in_mbr(p, query_mbr)) {
				out[cnt++] = p;
			}
		}
		return;
	}
	auto cur_inte = static_cast<inte_node *>(x.get());
	if (cur_inte->l_son != nullptr) {
		range_report_node(cur_inte->l_son, query_mbr, cnt, out);
	}
	if (cur_inte->r_son != nullptr) {
		range_report_node(cur_inte->r_son, query_mbr, cnt, out);
	}
}

template <class Out>
void Tree::range_report(bounding_box_type &query_mbr, size_t &cnt, Out &out)
{
	range_report_node(root, query_mbr, cnt, out);
}

auto Tree::knn_report(size_t &k, Point query_point)
{
	priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
	knn_report_node(root, k, query_point, nn_res);
	return nn_res;
}

template <class T>
void Tree::knn_report_node(shared_ptr<base_node> &x, size_t &k,
			   Point query_point, T &nn_res)
{
	if (!x)
		return;
	// print_mbr(x->mbr);
	if (x->is_leaf()) {
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		for (auto &p : cur_leaf->records) {
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
	auto cur_inte = static_cast<inte_node *>(x.get());
	auto l_son_sqrdis = ft_inf_max, r_son_sqrdis = ft_inf_max;

	if (cur_inte->l_son != nullptr) {
		l_son_sqrdis =
			point_mbr_sqrdis(query_point, cur_inte->l_son->mbr);
	}
	if (cur_inte->r_son != nullptr) {
		r_son_sqrdis =
			point_mbr_sqrdis(query_point, cur_inte->r_son->mbr);
	}

	if (l_son_sqrdis <= r_son_sqrdis) { // first go left
		if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second) {
			knn_report_node(cur_inte->l_son, k, query_point,
					nn_res);
		}
		if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second) {
			knn_report_node(cur_inte->r_son, k, query_point,
					nn_res);
		}
	} else { // first go right
		if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second) {
			knn_report_node(cur_inte->r_son, k, query_point,
					nn_res);
		}
		if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second) {
			knn_report_node(cur_inte->l_son, k, query_point,
					nn_res);
		}
	}
	return;
}

auto Tree::collect_records(shared_ptr<base_node> &x)
{
	sequence<Point> ret = {};
	if (!x) {
		return ret;
	}
	// cout << "[Cur MBR]: ";
	// print_mbr(x->mbr);
	if (x->is_leaf()) {
		// cout << "meet leaf node" << endl;
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		ret = cur_leaf->records;
		// for (auto &p: ret){
		//   cout << p << endl;
		// }
		return ret;
	}
	// cout << "meet inte node" << endl;
	auto cur_inte = static_cast<inte_node *>(x.get());
	// if (x == root){
	//   cout << "root sz: " << cur_inte->get_num_points() << endl;
	//   cout << "lson sz: " << cur_inte->l_son->get_num_points() << endl;
	//   cout << "rson sz: " << cur_inte->r_son->get_num_points() << endl;
	// }
	sequence<Point> R;
	auto collect_left = [&]() { ret = collect_records(cur_inte->l_son); };
	auto collect_right = [&]() { R = collect_records(cur_inte->r_son); };
	par_do_if(x->get_num_points() >= granularity_cutoff, collect_left,
		  collect_right);
	ret.append(R);
	return ret;
}

auto Tree::check_mbr(shared_ptr<base_node> &x)
{
	if (!x) {
		return empty_mbr;
	}
	if (x->is_leaf()) {
		auto cur_leaf = static_cast<leaf_node *>(x.get());
		auto mbr = get_mbr(cur_leaf->records);
		if (!is_same_mbr(mbr, x->mbr)) {
			cout << "[ERROR]: incorrect MBR" << endl;
		}
		return mbr;
	}

	auto cur_inte = static_cast<inte_node *>(x.get());
	geobase::bounding_box_type l_mbr, r_mbr;
	auto collect_left = [&]() { l_mbr = check_mbr(cur_inte->l_son); };
	auto collect_right = [&]() { r_mbr = check_mbr(cur_inte->r_son); };
	par_do_if(x->get_num_points() >= granularity_cutoff, collect_left,
		  collect_right);
	auto ret = merge_mbr(l_mbr, r_mbr);
	if (!is_same_mbr(ret, x->mbr)) {
		cout << "[ERROR]: incorrect MBR" << endl;
	}
	return ret;
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight = 6,
	  uint_fast8_t ImbaRatio = 30>
class zdtree
    : public psi::base_tree<Point,
			    zdtree<Point, SplitRule, SkHeight, ImbaRatio>,
			    SkHeight, ImbaRatio>
{
public:
	using base_type =
		psi::base_tree<Point,
			       zdtree<Point, SplitRule, SkHeight, ImbaRatio>,
			       SkHeight, ImbaRatio>;
	using bucket_type = typename base_type::bucket_type;
	using balls_type = typename base_type::balls_type;
	using dims_type = typename base_type::dims_type;
	using bucket_seq_type = typename base_type::bucket_seq_type;
	using ball_seq_type = typename base_type::ball_seq_type;
	using coord_type = typename Point::coord_type;
	using coords_type = typename Point::coords_type;
	using slice_type = typename base_type::slice_type;
	using points_type = typename base_type::points_type;
	using points_iter_type = typename base_type::points_iter_type;
	using box_type = typename base_type::box_type;
	using box_seq_type = typename base_type::box_seq_type;
	using circle_type = typename base_type::normal_circle;

	using hyper_plane_type = typename base_type::hyper_plane_type;
	using hyper_plane_seq_type = typename base_type::hyper_plane_seq_type;
	using node_tag_type = typename base_type::node_tag_type;
	using node_tag_seq_type = typename base_type::node_tag_seq_type;
	using node_box_type = typename base_type::node_box_type;
	using node_box_seq_type = typename base_type::node_box_seq_type;
	using leaf_type = psi::node;
	using interior_type = psi::node;

	void delete_tree()
	{
		tree.clear();
		return;
	}

	geobase::bounding_box_type largest_mbr = base_type::get_empty_box();
	// convert to zdtree point format, with storing Z-values
	template <typename Range>
	auto point_convert(Range &in)
	{
		slice_type A = parlay::make_slice(in);
		parlay::sequence<geobase::Point> P(in.size());
		// FT x_min(ft_inf_max), x_max(ft_inf_min), y_min(ft_inf_max),
		//     y_max(ft_inf_min);
		parlay::parallel_for(0, in.size(), [&](size_t i) {
			P[i] = geobase::Point(in[i].aug.id, in[i].pnt[0],
					      in[i].pnt[1], in[i].pnt[2]);
			// P[i].id = in[i].aug.id;
			// P[i].x = in[i].pnt[0];
			// P[i].y = in[i].pnt[1];
			// x_max = max(x_max, in[i].x);
			// x_min = min(x_min, in[i].x);
			// y_max = max(y_max, in[i].y);
			// y_min = min(y_min, in[i].y);
			// in[i].morton_id = in[i].interleave_bits();
			// P[i].morton_id = SplitRule::encode(in[i]); // order
			// is not xyzxyzxyz...
		});
		// largest_mbr = base_type::get_box(largest_mbr,
		// base_type::get_box(in)); largest_mbr =
		// geobase::bounding_box_type( {geobase::Point(x_min, y_min),
		// geobase::Point(x_max, y_max)});
		return in;
	}

	auto box_convert(box_type const &q)
	{
		// geobase::bounding_box_type cur_q(
		//     {geobase::Point(q.first.pnt[0], q.first.pnt[1]),
		//      geobase::Point(q.second.pnt[0], q.second.pnt[1])});
		// return cur_q;
		return q;
	}

	template <typename Pset>
	void check_first(Pset &P, size_t ed = 10)
	{
		std::cout << "sz: " << P.size() << std::endl;
		for (auto i = 0; i != ed; i++) {
			cout << P[i] << endl;
			// std::cout << fixed << setprecision(6) << P[i].id <<
			// ", " << P[i].x <<
			// ", "
			//           << P[i].y << ", " << P[i].z << ", " <<
			//           P[i].morton_id << std::endl;
		}
	}

	template <typename Range>
	bool check_pt(Range &in, geobase::bounding_box_type &q)
	{
		for (size_t i = 0; i < in.size(); i++) {
			if (geobase::point_in_mbr(in[i], q)) {
				return true;
			}
		}
		return false;
	}

	void run_tests()
	{
		struct test_case {
			bounding_box_type a, b;
			int expected;
			std::string description;
		};
		using geobase::Point;
		std::vector<test_case> tests = {
			// Basic cases
			{bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(3.0, 3.0, 3.0)}),
			 bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(5.0, 5.0, 5.0)}),
			 1, "Full Containment"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(4.0, 4.0, 4.0)}),
			 bounding_box_type({Point(2.0, 2.0, 2.0)},
					   {Point(6.0, 6.0, 6.0)}),
			 0, "Partial Overlap"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 bounding_box_type({Point(2.0, 0.0, 0.0)},
					   {Point(4.0, 2.0, 2.0)}),
			 0, "Touching at Face"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 bounding_box_type({Point(2.0, 2.0, 0.0)},
					   {Point(4.0, 4.0, 2.0)}),
			 0, "Touching at Edge"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(1.0, 1.0, 1.0)}),
			 bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 0, "Touching at Point"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(1.0, 1.0, 1.0)}),
			 bounding_box_type({Point(2.0, 0.0, 0.0)},
					   {Point(3.0, 1.0, 1.0)}),
			 -1, "Disjoint"},

			{bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(3.0, 3.0, 3.0)}),
			 bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(3.0, 3.0, 3.0)}),
			 1, "Equal Boxes"},

			// More complicated cases
			{bounding_box_type({Point(2.5, 3.5, 1.5)},
					   {Point(4.5, 5.5, 3.5)}),
			 bounding_box_type({Point(1.0, 2.0, 0.0)},
					   {Point(6.0, 7.0, 5.0)}),
			 1, "Containment not aligned with origin"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(3.0, 3.0, 1.0)}),
			 bounding_box_type({Point(1.0, 1.0, 2.0)},
					   {Point(4.0, 4.0, 3.0)}),
			 -1, "XY overlap but Z disjoint"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(3.0, 3.0, 1.0)}),
			 bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 0, "Thin slab overlap at z=1"},

			{bounding_box_type({Point(2.0, 2.0, 0.0)},
					   {Point(4.0, 4.0, 2.0)}),
			 bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(5.0, 5.0, 3.0)}),
			 0, "Nested overlap but not full containment"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(10.0, 10.0, 10.0)}),
			 bounding_box_type({Point(2.0, 2.0, 2.0)},
					   {Point(3.0, 3.0, 3.0)}),
			 0, "Large box passed as 'small_mbr'"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 bounding_box_type({Point(1.5, 1.5, 1.5)},
					   {Point(3.0, 3.0, 3.0)}),
			 0, "Diagonal partial overlap small cube"},

			{bounding_box_type({Point(0.0, 0.0, 0.0)},
					   {Point(1.0, 1.0, 1.0)}),
			 bounding_box_type({Point(1.0, 1.0, 1.0)},
					   {Point(2.0, 2.0, 2.0)}),
			 0, "Diagonal single point touch"},
		};

		for (size_t i = 0; i < tests.size(); i++) {
			int result = zy_mbr_relation(tests[i].a, tests[i].b);
			std::cout << "Case " << i + 1 << " ("
				  << tests[i].description << "): "
				  << "Expected = " << tests[i].expected
				  << ", Got = " << result
				  << (result == tests[i].expected ? " [OK]"
								  : " [FAIL]")
				  << "\n";
			assert(result == tests[i].expected);
		}

		std::cout << "✅ All tests passed!\n";
	}

	template <typename Range>
	void build(Range in)
	{
		// auto P = point_convert(in);
		// run_tests();
		auto p_set = geobase::get_sorted_points(in);
		// check_first(p_set, 12);
		tree.build(p_set);

		// cout << "check mbr start" << endl;
		// tree.check_mbr(tree.root);
		// cout << "check mbr end" << endl;

		// auto p_list = tree.collect_records(tree.root);
		// auto q =
		// geobase::bounding_box_type({geobase::Point(14.0, 2.0, 8.0),
		// geobase::Point(65536.0, 99992.0, 99997.0)}); size_t sz = 0;
		// auto Out =
		// parlay::sequence<geobase::Point>::uninitialized(12);
		// tree.range_report(q, largest_mbr, sz, Out, 3); cout << "sz =
		// " << sz << endl; cout << Out.size() << endl; for (auto &p:
		// Out){
		//   cout << p << endl;
		// }
		// // cout << "collected sz = " << p_list.size() << endl;
		// int xx; cin >> xx;
	}

	void batch_insert(slice_type in)
	{
		auto p_set = geobase::get_sorted_points(in);
		tree.batch_insert_sorted(p_set);

		// cout << "check mbr start" << endl;
		// tree.check_mbr(tree.root);
		// cout << "check mbr end" << endl;
	}

	void batch_delete(slice_type in)
	{
		auto p_set = geobase::get_sorted_points(in);
		// cout << "remove: " << in.size() << endl;
		// auto sz0 = tree.collect_records(tree.root).size();
		// cout << "before: " << sz0 << endl;
		tree.batch_delete_sorted(p_set);
		// auto sz1 = tree.collect_records(tree.root).size();
		// cout << "after: " << sz1 << endl;
		// cout << "check mbr start" << endl;
		// tree.check_mbr(tree.root);
		// cout << "check mbr end" << endl;
	}

	template <typename node, typename Range>
	auto knn(node *T, Point const &q, psi::bounded_queue<Point, Range> &bq)
	{
		// auto num_pts_in_tree =
		// tree.collect_records(tree.root).size(); cout <<
		// num_pts_in_tree << endl;
		// // print_mbr(tree.root->mbr);
		// if (num_pts_in_tree != 5000000){
		//   cout << "[ERROR] # of points in tree: " << num_pts_in_tree
		//   << endl;
		// }

		psi::knn_logger logger;
		// geobase::Point cnv_q(q[0], q[1]);
		size_t k = bq.max_size();

		auto knnsqrdis = tree.knn_report(k, q).top().second;

		// cout << knnsqrdis << endl;
		// auto all_pts = tree.collect_records(tree.root);
		// auto bfsqrdis = geobase::knn_bf(k, q, all_pts);
		// if (dcmp(knnsqrdis - bfsqrdis) != 0){
		//   cout << "[ERROR] Incorrect knn results" << endl;
		//   cout << bfsqrdis << ", " << knnsqrdis << endl;
		// }

		return logger;
	}

	auto range_count(box_type const &q)
	{
		// auto num_pts_in_tree =
		// tree.collect_records(tree.root).size(); if (num_pts_in_tree
		// != 5000000){
		//   cout << "[ERROR] # of points in tree: " << num_pts_in_tree
		//   << endl;
		// }

		psi::range_query_logger logger;
		int size = 0;
		auto cnv_q = box_convert(q);

		size = tree.range_count(cnv_q);

		return std::make_pair(size, logger);
	}

	template <typename Range>
	auto range_query(box_type const &q, Range &&Out)
	{
		// auto num_pts_in_tree =
		// tree.collect_records(tree.root).size(); if (num_pts_in_tree
		// != 5000000){
		//   cout << "[ERROR] # of points in tree: " << num_pts_in_tree
		//   << endl;
		// }

		psi::range_query_logger logger;
		auto cnv_q = box_convert(q);
		// std::cout << cnv_q.first.x << ", " << cnv_q.first.y << ", "
		//           << cnv_q.second.x << ", " << cnv_q.second.y <<
		//           std::endl;
		size_t size = 0;
		// parlay::sequence<geobase::Point> ret(Out.size());

		// auto p_list = tree.collect_records(tree.root);
		// auto check_q =
		// geobase::bounding_box_type({geobase::Point(65536.0, 2.0, 8.0),
		// geobase::Point(98532.0, 99992.0, 99997.0)}); std::cout <<
		// std::endl << "check pt: " << check_pt(p_list, check_q) <<
		// std::endl;

		tree.range_report(cnv_q, size, Out);

		// assert(Out.size() == size);
		// if (size != Out.size()){··
		//   std::cout << "[ERROR]: " << size << ", " << Out.size() <<
		//   std::endl;
		// }
		// std::cout << (Out.size() == size) << std::endl;
		return std::make_pair(size, logger);
	}

	constexpr static char const *get_tree_name()
	{
		return "ZdTree";
	}
	constexpr static char const *check_has_box()
	{
		return "WithoutBox";
	}
	Tree tree = Tree(32);
};

} // namespace ZD3D
