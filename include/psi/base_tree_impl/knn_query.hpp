#ifndef PSI_BASE_TREE_IMPL_KNN_QUERY_HPP
#define PSI_BASE_TREE_IMPL_KNN_QUERY_HPP

#include <algorithm>
#include <utility>

#include "psi/base_tree.h"
#include "parlay/primitives.h"
#include "psi/dependence/tree_node.h"

namespace psi
{

/*
 * distance between two points_type
 * TODO: change the name to p2p_distance_square to avoid ambiguous
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::dis_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::p2p_distance_square(
	Point const &p, Point const &q)
{
	constexpr uint_fast8_t num_dims = Point::get_dim();
	dis_type r = 0;

	if constexpr (num_dims == 2) {
		dis_type x = static_cast<dis_type>(p.pnt[0]) -
			     static_cast<dis_type>(q.pnt[0]);
		dis_type y = static_cast<dis_type>(p.pnt[1]) -
			     static_cast<dis_type>(q.pnt[1]);
		r = x * x + y * y;
	} else if constexpr (num_dims == 3) {
		dis_type x = static_cast<dis_type>(p.pnt[0]) -
			     static_cast<dis_type>(q.pnt[0]);
		dis_type y = static_cast<dis_type>(p.pnt[1]) -
			     static_cast<dis_type>(q.pnt[1]);
		dis_type z = static_cast<dis_type>(p.pnt[2]) -
			     static_cast<dis_type>(q.pnt[2]);
		r = x * x + y * y + z * z;
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i]));
		}
	}
	return r;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::dis_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::p2b_min_distance_square(
	Point const &p, typename base_tree<Point, DerivedTree, SkHeight,
					   ImbaRatio>::box_type const &a)
{
	dis_type r = 0;
	/* the distance is 0 when p is inside the box */
	for (dims_type i = 0; i < num_dims; ++i) {
		if (num_type::lt(p.pnt[i], a.first.pnt[i])) {
			r += (static_cast<dis_type>(a.first.pnt[i]) -
			      static_cast<dis_type>(p.pnt[i])) *
			     (static_cast<dis_type>(a.first.pnt[i]) -
			      static_cast<dis_type>(p.pnt[i]));
		} else if (num_type::gt(p.pnt[i], a.second.pnt[i])) {
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(a.second.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(a.second.pnt[i]));
		} else { /* will not count the dis if p is inside the box in */
			 /* dimension i */
			;
		}
	}
	return r;
}

/* max distance between a Point and a box_type */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::dis_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::p2b_max_distance_square(
	Point const &p, typename base_tree<Point, DerivedTree, SkHeight,
					   ImbaRatio>::box_type const &a)
{
	dis_type r = 0;
	for (dims_type i = 0; i < num_dims; ++i) {
		if (num_type::lt(p.pnt[i],
				 (a.second.pnt[i] + a.first.pnt[i]) / 2)) {
			r += (static_cast<dis_type>(a.second.pnt[i]) -
			      static_cast<dis_type>(p.pnt[i])) *
			     (static_cast<dis_type>(a.second.pnt[i]) -
			      static_cast<dis_type>(p.pnt[i]));
		} else {
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(a.first.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(a.first.pnt[i]));
		}
	}
	return r;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline double
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::p2c_min_distance(
	Point const &p, Point const &center, dis_type const r)
{

	return std::sqrt(p2p_distance_square(p, center)) -
	       static_cast<double>(r);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline double
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::p2c_min_distance(
	Point const &p, CircleType const &cl)
{
	return p2c_min_distance(p, cl.get_center(), cl.get_radius());
}

/*
 * early return the partial distance between p and q if it is larger than
 * r else return the distance between p and q
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::dis_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::interruptible_distance(
	Point const &p, Point const &q, dis_type up)
{
	dis_type r = 0;
	dims_type i = 0;
	if (num_dims >= 6) {
		while (1) {
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i]));
			++i;
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i]));
			++i;
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i]));
			++i;
			r += (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i])) *
			     (static_cast<dis_type>(p.pnt[i]) -
			      static_cast<dis_type>(q.pnt[i]));
			++i;

			if (num_comparator<dis_type>::gt(r, up)) {
				return r;
			}
			if (i + 4 > num_dims) {
				break;
			}
		}
	}
	while (i < num_dims) {
		r += (static_cast<dis_type>(p.pnt[i]) -
		      static_cast<dis_type>(q.pnt[i])) *
		     (static_cast<dis_type>(p.pnt[i]) -
		      static_cast<dis_type>(q.pnt[i]));
		++i;
	}
	return r;
}

/* knn search for Point q */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_leaf(
	node *T, Point const &q, bounded_queue<Point, Range> &bq)
{
	assert(T->is_leaf);

	leaf_type *tl = static_cast<leaf_type *>(T);
	size_t i = 0;
	while (!bq.full() && i < tl->size) {
		bq.insert(std::make_pair(
			std::ref(tl->pts[(!tl->is_dummy) * i]),
			p2p_distance_square(q, tl->pts[(!tl->is_dummy) * i])));
		i++;
	}
	while (i < tl->size) {
		auto r = interruptible_distance(q, tl->pts[(!tl->is_dummy) * i],
						bq.top_value());
		if (num_type::lt(r, bq.top_value())) { /* the queue is full, */
						       /*
							* no need to insert points
							* with equal distances
							*/
			bq.insert(std::make_pair(
				std::ref(tl->pts[(!tl->is_dummy) * i]), r));
		} else if (tl->is_dummy) {
			break;
		}
		i++;
	}
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_binary(
	node *T, Point const &q, bounded_queue<Point, Range> &bq,
	box_type const &node_box, knn_logger &logger)
	requires std::same_as<typename interior_type::st_type,
			      typename base_tree<Point, DerivedTree, SkHeight,
						 ImbaRatio>::hyper_plane_type>
{
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		knn_leaf<leaf_type>(T, q, bq);
		return;
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);
	bool go_left =
		num_type::gt(ti->split.first - q.pnt[ti->split.second], 0);
	box_cut_type box_cut(node_box, ti->split, go_left);
	logger.generate_box_num += 1;

	knn_binary<leaf_type, interior_type>(go_left ? ti->left : ti->right, q,
					     bq, box_cut.get_first_box_cut(),
					     logger);

	logger.check_box_num++;
	if (num_type::gt(
		    p2b_min_distance_square(q, box_cut.get_second_box_cut()),
		    bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	knn_binary<leaf_type, interior_type>(go_left ? ti->right : ti->left, q,
					     bq, box_cut.get_box(), logger);
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_binary_box(
	node *T, Point const &q, bounded_queue<Point, Range> &bq,
	knn_logger &logger)
{
	if (bq.size() &&
	    num_type::gt(p2b_min_distance_square(
				 q, retrieve_box<leaf_type, interior_type>(T)),
			 bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}

	if (T->is_leaf) {
		logger.vis_leaf_num++;
		knn_leaf<leaf_type>(T, q, bq);
		return;
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);
	coord_type dist_left = p2b_min_distance_square(
		q, retrieve_box<leaf_type, interior_type>(ti->left));
	coord_type dist_right = p2b_min_distance_square(
		q, retrieve_box<leaf_type, interior_type>(ti->right));
	bool go_left = num_type::leq(dist_left, dist_right);

	knn_binary_box<leaf_type, interior_type>(go_left ? ti->left : ti->right,
						 q, bq, logger);

	logger.check_box_num++;
	if (num_type::gt(go_left ? dist_right : dist_left, bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	knn_binary_box<leaf_type, interior_type>(go_left ? ti->right : ti->left,
						 q, bq, logger);
	return;
}

/* compute knn for multinode as if a binary node */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_multi_expand(
	node *T, Point const &q, dims_type dim, bucket_type idx,
	bounded_queue<Point, Range> &bq, box_type const &node_box,
	knn_logger &logger)
{
	if (T->size == 0) {
		return;
	}

	if (T->is_leaf) {
		logger.vis_leaf_num++;
		knn_leaf<leaf_type>(T, q, bq);
		return;
	}

	logger.vis_interior_num++;

	interior_type *ti = static_cast<interior_type *>(T);
	auto const &split =
		interior_type::equal_split() ? ti->split[dim] : ti->split[idx];
	bool go_left = num_type::gt(split.first - q.pnt[split.second], 0);

	bucket_type first_idx = (idx << 1) + static_cast<bucket_type>(!go_left);
	bucket_type second_idx = (idx << 1) + static_cast<bucket_type>(go_left);
	bool reach_leaf =
		first_idx >= interior_type::get_regions(); /* whether reach the
							      skeleton leaf */
	node *first_node =
		reach_leaf ? ti->tree_nodes[first_idx -
					    interior_type::get_regions()]
			   : T;
	node *second_node =
		reach_leaf ? ti->tree_nodes[second_idx -
					    interior_type::get_regions()]
			   : T;
	if (reach_leaf) {
		first_idx = second_idx = 1;
	}

	box_cut_type box_cut(node_box, split, go_left);
	logger.generate_box_num += 1;
	assert((dim + 1) % num_dims != 0 ||
	       (first_idx == 1 && second_idx == 1));

	knn_multi_expand<leaf_type, interior_type>(
		first_node, q, (dim + 1) % num_dims, first_idx, bq,
		box_cut.get_first_box_cut(), logger);

	/* compute the other bounding box */
	logger.check_box_num++;
	if (num_type::gt(
		    p2b_min_distance_square(q, box_cut.get_second_box_cut()),
		    bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	knn_multi_expand<leaf_type, interior_type>(
		second_node, q, (dim + 1) % num_dims, second_idx, bq,
		box_cut.get_box(), logger);
	return;
}

/* compute knn for multinode as if a binary node */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_multi_expand_box(
	node *T, Point const &q, dims_type dim, bucket_type idx,
	bounded_queue<Point, Range> &bq, knn_logger &logger)
{
	if (T->size == 0) {
		return;
	}

	if (idx == 1 && bq.size() &&
	    num_type::gt(p2b_min_distance_square(
				 q, retrieve_box<leaf_type, interior_type>(T)),
			 bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		knn_leaf<leaf_type>(T, q, bq);
		return;
	}

	logger.vis_interior_num++;

	interior_type *ti = static_cast<interior_type *>(T);
	auto const &split =
		interior_type::equal_split() ? ti->split[dim] : ti->split[idx];
	bool go_left = num_type::gt(split.first - q.pnt[split.second], 0);

	bucket_type first_idx = (idx << 1) + static_cast<bucket_type>(!go_left);
	bucket_type second_idx = (idx << 1) + static_cast<bucket_type>(go_left);
	bool reach_leaf =
		first_idx >= interior_type::get_regions(); /* whether reach the
							      skeleton leaf */
	node *first_node =
		reach_leaf ? ti->tree_nodes[first_idx -
					    interior_type::get_regions()]
			   : T;
	node *second_node =
		reach_leaf ? ti->tree_nodes[second_idx -
					    interior_type::get_regions()]
			   : T;
	if (reach_leaf) {
		first_idx = second_idx = 1;
	}
	auto second_box =
		second_node->is_leaf
			? retrieve_box<leaf_type, interior_type>(second_node)
			: static_cast<interior_type *>(second_node)
				  ->get_box_by_id(second_idx);

	assert((dim + 1) % num_dims != 0 ||
	       (first_idx == 1 && second_idx == 1));

	knn_multi_expand_box<leaf_type, interior_type>(
		first_node, q, (dim + 1) % num_dims, first_idx, bq, logger);

	/* compute the other bounding box */
	logger.check_box_num++;
	if (num_type::gt(p2b_min_distance_square(q, second_box),
			 bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	knn_multi_expand_box<leaf_type, interior_type>(
		second_node, q, (dim + 1) % num_dims, second_idx, bq, logger);
	return;
}

/* compute knn for multi-node by computing bounding boxes */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_multi(
	node *T, Point const &q, bounded_queue<Point, Range> &bq,
	knn_logger &logger)
{
	if (T->size == 0) {
		return;
	}

	if (bq.size() &&
	    num_type::gt(p2b_min_distance_square(
				 q, retrieve_box<leaf_type, interior_type>(T)),
			 bq.top_value()) &&
	    bq.full()) {
		logger.skip_box_num++;
		return;
	}
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		knn_leaf<leaf_type>(T, q, bq);
		return;
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);

	/* TODO: may change to use parlay::sequence */
	std::array<std::pair<coord_type, bucket_type>,
		   interior_type::get_regions()>
		dists;

	std::ranges::generate(dists, [i = 0, &q, ti]() mutable {
		auto r = std::make_pair(
			p2b_min_distance_square(
				q, retrieve_box<leaf_type, interior_type>(
					   ti->tree_nodes[i])),
			i);
		i++;
		return r;
	});
	std::ranges::sort(dists, std::less<>(),
			  [&](auto const &box_pair) { return box_pair.first; });

	knn_multi<leaf_type, interior_type>(ti->tree_nodes[dists[0].second], q,
					    bq, logger);
	for (bucket_type i = 1; i < interior_type::get_regions(); ++i) {
		logger.check_box_num++;
		if (num_type::gt(dists[i].first, bq.top_value()) && bq.full()) {
			logger.skip_box_num++;
			continue;
		}
		knn_multi<leaf_type, interior_type>(
			ti->tree_nodes[dists[i].second], q, bq, logger);
	}

	return;
}

/*
 * Convenience knn: sizes and fills the queue itself, then hands back copies
 * sorted by distance. The buffer-taking form is still there for callers that
 * want to own the storage.
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn_result_seq_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::knn(Point const &q,
							size_t k) const
{
	size_t const n =
		std::min(k, static_cast<DerivedTree const *>(this)->get_size());
	if (n == 0) {
		return knn_result_seq_type();
	}

	/* bounded_queue holds references, which have no default state, so the
	 * buffer needs a filler; it is overwritten before anything reads it. */
	using nn_pair = std::pair<std::reference_wrapper<Point>, dis_type>;
	Point filler = q;
	parlay::sequence<nn_pair> buf(n, nn_pair(std::ref(filler), dis_type()));
	bounded_queue<Point, nn_pair> bq(parlay::make_slice(buf));

	static_cast<DerivedTree const *>(this)->knn(q, bq);

	knn_result_seq_type out(bq.size());
	for (size_t i = 0; i < bq.size(); i++) {
		out[i] = knn_result_type(buf[i].first.get(), buf[i].second);
	}
	std::sort(out.begin(), out.end(),
		  [](knn_result_type const &a, knn_result_type const &b) {
			  return a.second < b.second;
		  });
	return out;
}

} /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_KNN_QUERY_HPP */
