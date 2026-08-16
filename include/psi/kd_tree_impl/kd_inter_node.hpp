#ifndef PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_
#define PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_

#include <cstddef>

#include "psi/kd_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi
{

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
struct kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	       ImbaRatio>::kd_interior_node
    : binary_node<Point, splitter_type, InteriorAugType> {
	using pt_type = Point;
	using st_type = splitter_type;
	using at_type = InteriorAugType;

	kd_interior_node(node *_left, node *_right, st_type const &_split)
	    : binary_node<Point, splitter_type, at_type>(
		      _left, _right, _split,
		      at_type(at_type::template create<leaf_type,
						       interior_type>(_left,
								      _right)))
	{
	}

	kd_interior_node(node *_left, node *_right, st_type const &_split,
			 at_type const &_aug)
	    : binary_node<Point, splitter_type, at_type>(_left, _right, _split,
							 _aug)
	{
	}

	// Adding a virtual destructor makes node polymorphic
	virtual ~kd_interior_node() = default;

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

	auto update_aug(node *l, node *r)
	{
		return this->aug.template update<leaf_type, interior_type>(l,
									   r);
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

	auto reset_aug()
	{
		return this->aug.reset();
	}
};

// template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
//           uint_fast8_t ImbaRatio>
// template <uint_fast8_t md>
// struct kd_tree<Point, SplitRule, SkHeight, ImbaRatio>::KdCompressionNode
//     : multi_node<Point, md, CompressNodeSplitter, AugType> {
//   using base_node = multi_node<Point, md, CompressNodeSplitter,
//   AugType>; using KdNodeArr = typename base_node::node_arr_type; using
//   pt_type = Point; using st_type = CompressNodeSplitter; using at_type =
//   AugType;

//   KdCompressionNode(KdNodeArr const& _tree_nodes, const st_type& _split,
//                     const at_type& _aug)
//       : base_node(_tree_nodes, _split, _aug) {}

//   // Adding a virtual destructor makes node polymorphic
//   virtual ~KdCompressionNode() = default;

//   inline void set_parallel_flag(bool const flag) { this->aug.emplace(flag); }

//   inline void reset_parallel_flag() { this->aug.reset(); }

//   inline bool get_parallel_flag_ini_status() { return this->aug.has_value();
//   }

//   // NOTE: use a tri-state bool to indicate whether a subtree needs to be
//   // rebuilt. If aug is not INITIALIZED, then it means there is no need to
//   // rebuild; otherwise, the value depends on the initial tree size before
//   // rebuilding.
//   inline bool force_parallel() const {
//     return this->aug.has_value() ? this->aug.value()
//                                  : this->size >
//                                  base_type::serial_build_cutoff;
//   }
// };

} // namespace psi

#endif // PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_
