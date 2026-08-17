#ifndef PSI_BASE_TREE_IMPL_BOX_CUT_HPP_
#define PSI_BASE_TREE_IMPL_BOX_CUT_HPP_

#include "psi/base_tree.h"

namespace psi
{
template <typename Traits, typename DerivedTree>
struct base_tree<Traits, DerivedTree>::box_cut_type {
	using base_type = base_tree<Traits, DerivedTree>;

	box_cut_type(box_type const &box, hyper_plane_type const &hp,
		     bool go_left)
	    : box(box), hp(hp), go_left(go_left)
	{
	}

	inline box_type const &get_first_box_cut()
	{
		mod_dim = go_left ? &box.second.pnt[hp.second]
				  : &box.first.pnt[hp.second];
		std::ranges::swap(hp.first, *mod_dim);
		return box;
	}

	inline box_type const &get_second_box_cut()
	{
		std::ranges::swap(hp.first, *mod_dim);
		mod_dim = go_left ? &box.first.pnt[hp.second]
				  : &box.second.pnt[hp.second];
		*mod_dim = hp.first;
		return box;
	}

	inline box_type const &get_box() const
	{
		return box;
	}

	inline hyper_plane_type const &get_hyper_plane() const
	{
		return hp;
	}

	box_type box;
	coord_type *mod_dim;
	hyper_plane_type hp; /* PARA: the split and the cutting dimension */
	bool const go_left;
};
}; /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_BOX_CUT_HPP_ */
