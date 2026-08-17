#ifndef PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_
#define PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_

#include <cstddef>

#include "psi/kd_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi
{

template <typename Traits>
struct kd_tree<Traits>::kd_interior_node
    : binary_node<point_type, splitter_type,
		  typename Traits::interior_aug_type> {
	using pt_type = point_type;
	using st_type = splitter_type;
	using at_type = typename Traits::interior_aug_type;

	kd_interior_node(node *_left, node *_right, st_type const &_split)
	    : binary_node<point_type, splitter_type, at_type>(
		      _left, _right, _split,
		      at_type(at_type::template create<leaf_type,
						       interior_type>(_left,
								      _right)))
	{
	}

	kd_interior_node(node *_left, node *_right, st_type const &_split,
			 at_type const &_aug)
	    : binary_node<point_type, splitter_type, at_type>(_left, _right,
							      _split, _aug)
	{
	}

	/* Adding a virtual destructor makes node polymorphic */
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

} /* namespace psi */

#endif /* PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_ */
