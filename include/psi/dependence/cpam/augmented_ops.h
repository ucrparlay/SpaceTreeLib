#pragma once

#include "map_ops.h"
#include "utils.h"

// *******************************************
//   AUGMENTED MAP OPERATIONS
// *******************************************

namespace psi
{
namespace cpam
{

template <class map_type>
struct augmented_ops : map_type {
	using Entry = typename map_type::Entry;
	using node = typename map_type::node;
	using Seq = typename map_type::_Seq;
	using et_type = typename map_type::et_type;
	using gc_type = typename map_type::gc_type;
	using K = typename map_type::K;
	using aug_t = typename Entry::aug_t;
	using ptr = typename gc_type::ptr;
	using map_type::B;
	using map_type::base_case_size;
	using map_type::node_limit;

	static inline aug_t aug_val(node *b)
	{
		return Seq::aug_val(b);
	}

	static inline aug_t const &aug_val_ref(node *b)
	{
		return Seq::aug_val_ref(b);
	}

	struct aug_sum_t {
		aug_t result;
		aug_sum_t() : result(Entry::get_empty())
		{
		}
		void add_entry(et_type e)
		{
			result = Entry::combine(result, Entry::from_entry(e));
		}
		void add_aug_val(aug_t av)
		{
			result = Entry::combine(result, av);
		}
	};

	// the sum left of or at key
	template <class aug>
	static void aug_sum_left(node *b, K const &key, aug &a)
	{
		while (b) {
			if (map_type::is_compressed(b)) {
				auto fn = [&](auto const &et) {
					if (map_type::comp(Entry::get_key(et),
							   key)) {
						a.add_entry(et);
						return true;
					}
					return false;
				};
				map_type::iterate_cond(b, fn);
				return;
			}
			auto rb = map_type::cast_to_regular(b);
			if (!map_type::comp(key, map_type::get_key(rb))) {
				a.add_entry(map_type::get_entry(rb));
				if (rb->lc)
					a.add_aug_val(
						map_type::aug_val(rb->lc));
				b = rb->rc;
			} else {
				b = rb->lc;
			}
		}
	}

	// the sum right of or at key
	template <class aug>
	static void aug_sum_right(node *b, K const &key, aug &a)
	{
		while (b) {
			if (map_type::is_compressed(b)) {
				auto fn = [&](auto const &et) {
					if (!map_type::comp(Entry::get_key(et),
							    key)) {
						a.add_entry(et);
					}
				};
				map_type::iterate_seq(b, fn);
				return;
			}
			auto rb = map_type::cast_to_regular(b);
			if (!map_type::comp(map_type::get_key(rb), key)) {
				a.add_entry(map_type::get_entry(rb));
				if (rb->rc)
					a.add_aug_val(
						map_type::aug_val(rb->rc));
				b = rb->lc;
			} else
				b = rb->rc;
		}
	}

	template <class aug>
	static void aug_sum_range(node *b, K const &key_left,
				  K const &key_right, aug &a)
	{
		node *r = map_type::range_root_2(b, key_left, key_right);
		if (r) {
			if (map_type::is_compressed(r)) {
				auto fn = [&](auto const &et) {
					if (map_type::comp(
						    key_left,
						    Entry::get_key(et)) &&
					    map_type::comp(Entry::get_key(et),
							   key_right)) {
						a.add_entry(et);
					}
				};
				map_type::iterate_seq(r, fn);
			} else {
				auto rr = map_type::cast_to_regular(r);
				// add in left side (right of or at key_left)
				aug_sum_right(rr->lc, key_left, a);
				// add in middle
				a.add_entry(map_type::get_entry(rr));
				// add in right side (left of or at key_right)
				aug_sum_left(rr->rc, key_right, a);
			}
		}
	}

	template <typename Func>
	static std::optional<et_type> aug_select(node *b, Func const &f)
	{
		if (!b)
			return {};
		if (map_type::is_compressed(b)) {
			std::optional<et_type> ret;
			auto fn = [&](auto const &et) {
				if (!f(Entry::from_entry(et))) {
					ret = et;
					return false; // stop
				}
				return true; // keep iterating
			};
			map_type::iterate_cond(b, fn);
			return ret;
		}
		auto rb = map_type::cast_to_regular(b);
		if (f(map_type::aug_val(rb->lc))) {
			if (f(Entry::from_entry(map_type::get_entry(rb)))) {
				return aug_select(rb->rc, f);
			}
			return map_type::get_entry(rb);
		}
		return aug_select(rb->lc, f);
	}

	template <class Func>
	static node *aug_filter_bc(ptr b1, Func const &f)
	{
		assert(b1.size() > 0);
		et_type stack[base_case_size + 1];

		auto b1_node = b1.node_ptr();
		size_t offset = 0;
		aug_t cur = Entry::get_empty();
		auto copy_f = [&](et_type a) { // has to be a copy since we move
			cur = Entry::combine(cur, Entry::from_entry(a));
			if (f(cur)) {
				parlay::move_uninitialized(stack[offset++], a);
			}
		};
		map_type::iterate_seq(b1_node, copy_f);
		assert(offset <= base_case_size);

		map_type::decrement_recursive(b1_node);

		if (offset < B) {
			return map_type::to_tree_impl((et_type *)stack, offset);
		} else {
			return map_type::make_compressed(stack, offset);
		}
	}

	template <class Func>
	static node *aug_filter(ptr b, Func const &f,
				size_t granularity = node_limit)
	{
		if (b.empty())
			return NULL;
		if (b.size() <= base_case_size) {
			return aug_filter_bc(std::move(b), f);
		}
		// TODO: better functionality for getting aug_val from b
		// std::cout << "My aug_val = " << aug_val(b.unsafe_ptr()) <<
		// std::endl;
		if (!f(aug_val(b.unsafe_ptr())))
			return NULL;

		size_t n = b.size();
		auto [lc, e, rc, root] = map_type::expose(std::move(b));

		auto [l, r] = utils::fork<node *>(
			n >= granularity,
			[&]() {
				return aug_filter(std::move(lc), f,
						  granularity);
			},
			[&]() {
				return aug_filter(std::move(rc), f,
						  granularity);
			});

		if (f(Entry::from_entry(e))) {
			return map_type::join(l, e, r, root);
		} else {
			gc_type::decrement(root);
			return map_type::join2(l, r);
		}
	}

	template <class Func>
	static node *insert_lazy(node *b, et_type const &e, Func const &f)
	{
		aug_t av = Entry::from_entry(e);
		auto g = [&](aug_t const &a) { return Entry::combine(av, a); };

		auto lazy_join = [&](node *l, node *r, node *_m) -> node * {
			auto m = map_type::cast_to_regular(_m);
			m->rc = r;
			m->lc = l;
			if (map_type::is_balanced(m)) {
				map_type::lazy_update(m, g);
				return m;
			} else
				return map_type::node_join(l, r, m);
		};

		return map_type::template insert_tmpl<Func, decltype(lazy_join),
						      false>(b, e, f,
							     lazy_join);
	}

	// F for box check, F2 for point check
	// template <typename F, typename F2>
	// static size_t range_count_filter2(node* b, F const& f, const F2& f2,
	//                                   size_t granularity = node_limit) {
	//   if (!b) return 0;
	//   auto cur_aug = aug_val(b);
	//   auto flag = f(cur_aug.first);
	//   if (flag < 0) return 0;  // exclude
	//   if (flag == 1) {
	//     return cur_aug.second;  // fully contained
	//   }
	//
	//   if (map_type::is_compressed(b)) {  // leaf node
	//     auto ret = 0;
	//     auto f_filter = [&](auto const& et) {
	//       auto cur_pt = std::get<1>(et);
	//       if (f2(cur_pt) == 1) {
	//         ret++;
	//       }
	//     };
	//     map_type::iterate_seq(b, f_filter);
	//     return ret;
	//   }
	//
	//   auto rb = map_type::cast_to_regular(b);
	//   auto cur_pt = map_type::get_val(rb);
	//   auto flag2 = f2(cur_pt) == 1 ? 1 : 0;
	//
	//   auto l = range_count_filter2(rb->lc, f, f2, granularity);
	//   auto r = range_count_filter2(rb->rc, f, f2, granularity);
	//
	//   return l + r + flag2;
	// }

	template <class base_tree, typename box_type, typename Logger>
	static size_t range_count_filter2(node *b, box_type const &query_box,
					  Logger &logger)
	{
		using base_type = base_tree;

		if (!b)
			return 0;

		auto const &node_box = map_type::aug_val_ref(b);
		if (!base_type::box_intersect_box(node_box, query_box)) {
			logger.skip_box_num++;
			return 0;
		} else if (base_type::within_box(node_box, query_box)) {
			logger.full_box_num++;
			return map_type::size(b); // fully contained
		}

		if (map_type::is_compressed(b)) { // leaf node
			logger.vis_leaf_num++;
			size_t ret = 0;
			auto f_filter = [&](auto const &et) {
				if (base_type::within_box(et, query_box)) {
					ret++;
				}
			};
			map_type::template iterate_seq<decltype(f_filter),
						       false>(b, f_filter);
			return ret;
		}

		logger.vis_interior_num++;
		auto rb = map_type::cast_to_regular(b);

		auto l = range_count_filter2<base_tree>(rb->lc, query_box,
							logger);
		auto r = range_count_filter2<base_tree>(rb->rc, query_box,
							logger);

		return l + r +
		       static_cast<size_t>(
			       base_type::within_box(rb->entry.first, query_box)
				       ? 1
				       : 0);
	}

	//  F is point-point dis, F2 is point-mbr dis
	template <class base_tree, typename F, typename F2, typename Out,
		  typename Logger>
	static void knn_filter(node *b, et_type const q, F const &f,
			       const F2 &f2, size_t &k, Out &out,
			       Logger &logger)
	{
		using base_type = base_tree;
		if (!b)
			return;

		auto pt_check = [&](auto &cur_pt) {
			auto cur_dis = f(cur_pt);
			// if (out.size() < k)
			if (!out.full())
				out.insert(std::make_pair(std::ref(cur_pt),
							  cur_dis));
			else if (cur_dis < out.top_value()) {
				// out.pop();
				out.insert(std::make_pair(std::ref(cur_pt),
							  cur_dis));
			}
		};

		if (map_type::is_compressed(b)) { // leaf node·
			logger.vis_leaf_num++;
			auto f_filter = [&](et_type &et) { pt_check(et); };
			map_type::iterate_seq(b, f_filter);
			return;
		}

		logger.vis_interior_num++;
		auto rb = map_type::cast_to_regular(b);
		// auto cur_pt = map_type::get_val(rb);
		auto cur_pt = rb->entry.first;
		pt_check(cur_pt);

		auto l_dis = std::numeric_limits<long>::max();
		auto r_dis = std::numeric_limits<long>::max();
		if (rb->lc) {
			auto cur_aug = aug_val(rb->lc);
			l_dis = f2(cur_aug);
		}
		if (rb->rc) {
			auto cur_aug = aug_val(rb->rc);
			r_dis = f2(cur_aug);
		}
		auto go_left = [&]() {
			if (!out.full() || l_dis < out.top().second) {
				knn_filter<base_type>(rb->lc, q, f, f2, k, out,
						      logger);
			}
		};
		auto go_right = [&]() {
			if (!out.full() || r_dis < out.top().second) {
				knn_filter<base_type>(rb->rc, q, f, f2, k, out,
						      logger);
			}
		};

		if (l_dis <= r_dis) { //  go left first
			go_left();
			go_right();
		} else {
			go_right();
			go_left();
		}
	}

	template <class base_tree, typename Logger, typename bounded_queue>
	static void knn(node *b, et_type const &q, bounded_queue &bq,
			Logger &logger)
	{
		using base_type = base_tree;
		using coord_type = typename et_type::coord_type;
		using dis_type = typename et_type::dis_type;

		if (!b)
			return;

		if (bq.size() &&
		    base_type::p2b_min_distance_square(
			    q, map_type::aug_val_ref(b)) > bq.top_value() &&
		    bq.full()) {
			logger.skip_box_num++;
			return;
		}

		if (map_type::is_compressed(b)) { // leaf node
			logger.vis_leaf_num++;
			auto f_filter = [&](et_type &et) {
				if (!bq.full()) {
					bq.insert(std::make_pair(
						std::ref(et),
						base_type::p2p_distance_square(
							q, et)));
					return;
				}
				auto r = base_type::interruptible_distance(
					q, et, bq.top_value());
				if (r <
				    bq.top_value()) { // PERF: remember
						      // currently the queue is
						      // full; if r == bq.top(),
						      // then it is useless to
						      // insert it, as it should
						      // not appears in the
						      // queue
					bq.insert(std::make_pair(std::ref(et),
								 r));
				}
				return;
			};
			map_type::template iterate_seq<decltype(f_filter),
						       false>(b, f_filter);
			return;
		}

		logger.vis_interior_num++;
		auto rb = map_type::cast_to_regular(b);
		dis_type d_lc =
			rb->lc ? base_type::p2b_min_distance_square(
					 q, map_type::aug_val_ref(rb->lc))
			       : std::numeric_limits<dis_type>::max();
		dis_type d_rc =
			rb->rc ? base_type::p2b_min_distance_square(
					 q, map_type::aug_val_ref(rb->rc))
			       : std::numeric_limits<dis_type>::max();
		bool go_left = d_lc <= d_rc;

		// check current entry
		// the rb->entry is a <point, aug> pair
		auto r = base_type::interruptible_distance(q, rb->entry.first,
							   bq.top_value());
		if (!bq.full() || r < bq.top_value()) {
			bq.insert(std::make_pair(std::ref(rb->entry.first), r));
		}

		knn<base_type>(go_left ? rb->lc : rb->rc, q, bq, logger);

		logger.check_box_num++;
		if (((go_left ? d_rc : d_lc) > bq.top_value()) && bq.full()) {
			logger.skip_box_num++;
			return;
		}

		knn<base_type>(go_left ? rb->rc : rb->lc, q, bq, logger);

		return;
	}

	template <class F, typename F2>
	static size_t range_count_filter(ptr b, F const &f, const F2 &f2,
					 size_t granularity = node_limit)
	{
		if (b.empty())
			return 0;
		// auto cur_par = map_type::get_entry(b.unsafe_ptr());
		// map_type::print_node_info(b.unsafe_ptr(), "cur");
		// std::cout << cur_val.first << "-" << cur_val.second <<
		// std::endl; auto cur_pt =
		// std::get<1>(map_type::get_entry(b.unsafe_ptr()));
		auto cur_aug = aug_val(b.unsafe_ptr());
		auto [lc, e, rc, root] = map_type::expose(std::move(b));
		// auto [lc2, e, rc, root] = map_type::expose(std::move(b));

		auto cur_pt = std::get<1>(e);
		// std::cout << std::fixed << std::setprecision(6) << cur_pt.x
		// << ", " << cur_pt.y << std::endl;

		auto flag = f(cur_aug.first);

		if (flag < 0) {
			gc_type::decrement(root);
			return 0;
		}
		if (flag == 1) {
			gc_type::decrement(root);
			// std::cout << "found " << cur_aug.second << " points"
			// << std::endl;
			return cur_aug.second;
		}

		// auto pt_box = std::make_pair(cur_pt, cur_pt);
		// auto cur_pt_inside = f(pt_box, 0) > 0 ? 1 : 0;
		auto cur_pt_inside = f2(cur_pt) > 0 ? 1 : 0;

		// size_t n = b.size();
		// auto [lc, e, rc, root] = map_type::expose(std::move(b));

		// auto [l, r] = utils::fork<size_t>(n >= granularity,
		//   [&]() {return range_count_filter(std::move(lc), f,
		//   granularity);},
		//   [&]() {return range_count_filter(std::move(rc), f,
		//   granularity);});

		auto l = range_count_filter(std::move(lc), f, f2, granularity);
		auto r = range_count_filter(std::move(rc), f, f2, granularity);

		gc_type::decrement(root);

		return l + r + cur_pt_inside;
	}

	// template <class F, typename Out>
	// static void range_report_filter2(node* b, F const& f, int64_t& cnt,
	// Out& out,
	//                                  size_t granularity = node_limit) {
	//   if (!b) return;
	//
	//   auto cur_aug = aug_val(b);
	//   auto flag = f(cur_aug.first);
	//   if (flag < 0) return;  // exclude
	//
	//   if (map_type::is_compressed(b)) {  // leaf node
	//     if (flag == 1) {
	//       auto f_filter = [&](auto const& et) { out[cnt++] =
	//       std::get<1>(et);
	//       }; map_type::iterate_seq(b, f_filter); return;  // fully
	//       contained
	//     }
	//
	//     auto f_filter = [&](auto const& et) {
	//       auto cur_pt = std::get<1>(et);
	//       auto pt_box = std::make_pair(cur_pt, cur_pt);
	//       if (f(pt_box) == 1) {
	//         out[cnt++] = cur_pt;
	//       }
	//     };
	//     map_type::iterate_seq(b, f_filter);
	//     return;
	//   }
	//
	//   auto rb = map_type::cast_to_regular(b);
	//   auto cur_pt = map_type::get_val(rb);
	//   auto pt_box = std::make_pair(cur_pt, cur_pt);
	//   auto flag2 = f(pt_box) == 1 ? 1 : 0;
	//   if (flag2) out[cnt++] = cur_pt;
	//
	//   range_report_filter2(rb->lc, f, cnt, out, granularity);
	//   range_report_filter2(rb->rc, f, cnt, out, granularity);
	// }

	// WARN: the output is not necessary sorted
	template <typename R>
	static size_t flatten(node *b, R out)
	{
		if (!b)
			return 0;
		if (map_type::is_compressed(b)) {
			size_t sz = 0;
			auto copy_f = [&](auto const &et) { out[sz++] = et; };
			map_type::template iterate_seq<decltype(copy_f), false>(
				b, copy_f);
			return sz;
		}
		auto rb = map_type::cast_to_regular(b);
		auto ls = map_type::size(rb->lc);
		assert(out.size() == rb->s);
		assert(rb->s ==
		       map_type::size(rb->lc) + map_type::size(rb->rc) + 1);

		utils::fork_no_result(
			rb->s >= 1024,
			[&]() { flatten(rb->lc, out.cut(0, ls)); },
			[&]() {
				flatten(rb->rc, out.cut(ls + 1, out.size()));
			});
		out[ls] = rb->entry.first;

		return rb->s;
	}

	template <class base_tree, typename box_type, typename Logger,
		  typename Out>
	static void range_report_filter2(node *b, box_type const &query_box,
					 size_t &cnt, Out out, Logger &logger)
	{
		using base_type = base_tree;

		if (!b)
			return;

		auto const &node_box = map_type::aug_val_ref(b);
		if (!base_type::box_intersect_box(node_box, query_box)) {
			logger.skip_box_num++;
			return;
		} else if (base_type::within_box(node_box, query_box)) {
			logger.full_box_num++;
			cnt += flatten(b,
				       out.cut(cnt, cnt + map_type::size(b)));
			return;
		}

		if (map_type::is_compressed(b)) { // leaf node
			logger.vis_leaf_num++;
			auto f_filter = [&](auto const &et) {
				if (base_type::within_box(et, query_box)) {
					out[cnt++] = et;
				}
			};
			map_type::template iterate_seq<decltype(f_filter),
						       false>(b, f_filter);
			return;
		}

		logger.vis_interior_num++;
		auto rb = map_type::cast_to_regular(b);

		if (base_type::within_box(rb->entry.first, query_box)) {
			out[cnt++] = rb->entry.first;
		}

		range_report_filter2<base_tree>(rb->lc, query_box, cnt, out,
						logger);
		range_report_filter2<base_tree>(rb->rc, query_box, cnt, out,
						logger);
	}

	template <class F, typename Out>
	static void range_report_filter(ptr b, F const &f, int64_t &cnt,
					Out &out,
					size_t granularity = node_limit)
	{
		if (b.empty())
			return;
		auto cur_aug = aug_val(b.unsafe_ptr());

		auto [lc, e, rc, root] = map_type::expose(std::move(b));

		auto flag = f(cur_aug.first);

		if (flag < 0) {
			gc_type::decrement(root);
			return; //  exclude
		}

		auto cur_pt = std::get<1>(e);
		auto pt_box = std::make_pair(cur_pt, cur_pt);
		auto cur_pt_inside = f(pt_box) > 0 ? 1 : 0;

		range_report_filter(std::move(lc), f, cnt, out);
		if (cur_pt_inside) {
			out[cnt++] = cur_pt;
		}
		range_report_filter(std::move(rc), f, cnt, out);

		gc_type::decrement(root);
		return;
	}
};

} // namespace cpam
} // namespace psi
