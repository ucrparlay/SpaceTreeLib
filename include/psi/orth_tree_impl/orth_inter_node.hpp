#ifndef PSI_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_

#include "../orth_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
struct orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		 ImbaRatio>::orth_interior_node
    : multi_node<Point, md, splitter_type, InteriorAugType> {
	using base_node = multi_node<Point, md, splitter_type, InteriorAugType>;
	using orth_node_arr_type = typename base_node::node_arr_type;
	using pt_type = Point;
	using st_type = splitter_type;
	using at_type = InteriorAugType;

	orth_interior_node(orth_node_arr_type const &_tree_nodes,
			   st_type const &_split)
	    : base_node(_tree_nodes, _split,
			at_type(at_type::template create<
				leaf_type, interior_type>(_tree_nodes)))
	{
	}

	orth_interior_node(orth_node_arr_type const &_tree_nodes,
			   st_type const &_split, at_type const &_aug)
	    : base_node(_tree_nodes, _split, _aug)
	{
	}

	constexpr auto get_sub_tree_num() const
	{
		return base_node::get_regions();
	}

	inline void set_parallel_flag(bool const flag)
	{
		this->aug.set_parallel_flag(flag);
	}

	inline void reset_parallel_flag()
	{
		this->aug.reset_parallel_flag();
	}

	inline bool get_parallel_flag_ini_status()
	{
		return this->aug.get_parallel_flag_ini_status();
	}

	inline bool force_parallel() const
	{
		return this->aug.force_parallel(this->size);
	}

	auto update_aug(orth_node_arr_type const &tree_nodes)
	{
		return this->aug.template update<leaf_type, interior_type>(
			tree_nodes);
	}

	decltype(auto) get_box()
		requires has_box<at_type>
	{
		return this->aug.get_box();
	}

	decltype(auto) get_box() const
		requires has_box<at_type>
	{
		return this->aug.get_box();
	}

	auto get_box_by_id_recursive(int idx)
		requires has_box<at_type>
	{
		if (idx >= base_node::get_regions()) {
			auto o = this->tree_nodes[idx -
						  base_node::get_regions()];
			return base_type::template retrieve_box<leaf_type,
								interior_type>(
				o);
		}
		return base_type::get_box(get_box_by_id_recursive(2 * idx),
					  get_box_by_id_recursive(2 * idx + 1));
	}

	auto get_box_by_id(int idx)
		requires has_box<at_type>
	{
		if (idx == 1) {
			return this->aug.get_box();
		}
		return get_box_by_id_recursive(idx);
	}

	auto reset_aug()
	{
		return this->aug.reset();
	}

	// NOTE: specific part
	// NOTE: compute the spliiter
	template <typename box_type>
	static st_type compute_splitter(box_type const &box)
	{
		st_type split;

		if constexpr (md == 2) { // for dim = 2
			split[0] = hyper_plane_type(
				base_type::get_box_mid(0, box), 0);
			split[1] = hyper_plane_type(
				base_type::get_box_mid(1, box), 1);
		} else if constexpr (md == 3) { // for dim = 3
			split[0] = hyper_plane_type(
				base_type::get_box_mid(0, box), 0);
			split[1] = hyper_plane_type(
				base_type::get_box_mid(1, box), 1);
			split[2] = hyper_plane_type(
				base_type::get_box_mid(2, box), 2);
		} else {
			for (dims_type i = 0; i < md; ++i) { // for dim > 3
				split[i] = hyper_plane_type(
					base_type::get_box_mid(i, box), i);
			}
		}

		return std::move(split);
	}

	// NOTE: Generate the hyperplane sequences for the node
	template <typename hyper_plane_seq_type>
	void generate_hyper_plane_seq(hyper_plane_seq_type &hyper_seq, auto idx,
				      auto deep)
	{
		if (idx >= base_node::get_regions()) {
			return;
		}
		hyper_seq[idx] = this->split[deep];
		generate_hyper_plane_seq(hyper_seq, idx * 2, deep + 1);
		generate_hyper_plane_seq(hyper_seq, idx * 2 + 1, deep + 1);
		return;
	}

	// NOTE: given an idx in the tree skeleton, modify its value its
	// position in left or right of the splitter
	static bucket_type sieve_point(Point const &p, st_type const &split,
				       bucket_type idx)
	{
		if constexpr (md == 2) { // for dim = 2
			idx = 2 * idx + 1 -
			      static_cast<bucket_type>(num_type::lt(
				      p.pnt[split[0].second], split[0].first));
			idx = 2 * idx + 1 -
			      static_cast<bucket_type>(num_type::lt(
				      p.pnt[split[1].second], split[1].first));
		} else if constexpr (md == 3) { // for dim = 3
			idx = 2 * idx + 1 -
			      static_cast<bucket_type>(num_type::lt(
				      p.pnt[split[0].second], split[0].first));
			idx = 2 * idx + 1 -
			      static_cast<bucket_type>(num_type::lt(
				      p.pnt[split[1].second], split[1].first));
			idx = 2 * idx + 1 -
			      static_cast<bucket_type>(num_type::lt(
				      p.pnt[split[2].second], split[2].first));
		} else {
			for (bucket_type i = 0; i < md; ++i) { // for dim > 3
				idx = 2 * idx + 1 -
				      static_cast<bucket_type>(num_type::lt(
					      p.pnt[split[i].second],
					      split[i].first));
			}
		}
		return idx;
	}

	bucket_type sieve_point(Point const &p, bucket_type idx)
	{
		return sieve_point(p, this->split, idx);
	}

	// NOTE: Given a box, produce new sub-boxes equivalent to modify the
	// existing boxes
	template <typename BoxSeqSlice>
	static void compute_subregions_rec(BoxSeqSlice box_seq,
					   st_type const &split,
					   bucket_type deep)
	{
		if (deep == md) {
			assert(box_seq.size() == 1);
			return;
		}

		bucket_type n = box_seq.size();
		assert(n - n / 2 == n / 2);

		std::for_each_n(box_seq.begin(), n / 2, [&](auto &bx) {
			bx.second.pnt[split[deep].second] = split[deep].first;
		});
		std::for_each_n(box_seq.begin() + n / 2, n - n / 2,
				[&](auto &bx) {
					bx.first.pnt[split[deep].second] =
						split[deep].first;
				});

		compute_subregions_rec(box_seq.cut(0, n / 2), split, deep + 1);
		compute_subregions_rec(box_seq.cut(n / 2, n), split, deep + 1);
		return;
	}

	template <typename box_seq_type, typename box_type>
	box_seq_type compute_subregions(box_type const &box)
	{
		auto box_seq = box_seq_type(base_node::get_regions(), box);
		compute_subregions_rec(parlay::make_slice(box_seq), this->split,
				       0);
		return std::move(box_seq);
	}

	template <typename box_seq_type, typename box_type>
	static box_seq_type compute_subregions(box_type const &box,
					       st_type const &split)
	{
		auto box_seq = box_seq_type(base_node::get_regions(), box);
		compute_subregions_rec(parlay::make_slice(box_seq), split, 0);
		return std::move(box_seq);
	}

	// NOTE: Given a box and a bucket id, modify the box for that bucket
	template <typename box_type>
	static void modify_box_by_id(bucket_type id, st_type const &split,
				     box_type &box)
	{
		assert(id >= 0 && id <= base_node::get_regions());

		if constexpr (md == 2) { // dim =2
			auto &target1 =
				(id & (1 << 1)) ? box.first : box.second;
			target1.pnt[split[0].second] = split[0].first;

			auto &target2 =
				(id & (1 << 0)) ? box.first : box.second;
			target2.pnt[split[1].second] = split[1].first;
		} else if constexpr (md == 3) { // dim =3
			auto &target1 =
				(id & (1 << 2)) ? box.first : box.second;
			target1.pnt[split[0].second] = split[0].first;

			auto &target2 =
				(id & (1 << 1)) ? box.first : box.second;
			target2.pnt[split[1].second] = split[1].first;

			auto &target3 =
				(id & (1 << 0)) ? box.first : box.second;
			target3.pnt[split[2].second] = split[2].first;
		} else { // dim > 3
			for (bucket_type i = md; i > 0; --i) {
				auto &target = (id & (1 << (i - 1)))
						       ? box.first
						       : box.second;
				target.pnt[split[md - i].second] =
					split[md - i].first;
			}
		}
		return;
	}

	template <typename box_type>
	void modify_box_by_id(bucket_type id, box_type &box)
	{
		modify_box_by_id(id, this->split, box);
		return;
	}

	// NOTE: return a new box by modifying the box for a specific region
	template <typename box_type>
	static box_type get_box_by_region_id(bucket_type id,
					     st_type const &split,
					     box_type const &box)
	{
		box_type new_box(box);
		modify_box_by_id(id, split, new_box);
		return std::move(new_box);
	}

	template <typename box_type>
	box_type get_box_by_region_id(bucket_type id, box_type const &box)
	{
		box_type new_box(box);
		modify_box_by_id(id, this->split, new_box);
		return std::move(new_box);
	}
};

} // namespace psi

#endif // PSI_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_
