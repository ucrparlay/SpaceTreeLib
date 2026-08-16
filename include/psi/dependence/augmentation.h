#ifndef PSI_DEPENDENCE_AUGMENTATION_H_
#define PSI_DEPENDENCE_AUGMENTATION_H_

#include <cstddef>
#include <optional>

#include "psi/dependence/tree_node.h"

namespace psi
{

/*
 * Default augmentations. kd_tree and orth_tree use these unless told otherwise,
 * so the common case needs no user code.
 *
 * base_type is a base_tree used only for its box vocabulary, so
 * base_tree<Point> works even though the tree instantiates base_tree<Point,
 * TheTree, ...>. Do not reach for anything here that varies with the other
 * parameters: bucket_type changes at SkHeight > 7.
 */

/* Bounding box of the points in a leaf. */
template <typename base_type>
struct box_leaf_aug {
	using box_type = typename base_type::box_type;
	using slice_type = typename base_type::slice_type;

	box_leaf_aug() : box(base_type::get_empty_box())
	{
	}
	box_leaf_aug(slice_type in) : box(base_type::get_box(in))
	{
	}

	box_type const &get_box() const
	{
		return box;
	}

	void update_aug(slice_type in)
	{
		box = base_type::get_box(in);
	}

	void reset()
	{
		box = base_type::get_empty_box();
	}

	box_type box;
};

/*
 * Bounding box of a subtree, plus the flag that picks parallel granularity.
 *
 * create and update are member templates because they are handed the very
 * interior_type that holds this augmentation. Both a binary and a
 * multi-way form are present; the node calls create<leaf_type,
 * interior_type>(...) and the last parameter is deduced from the argument.
 */
template <typename base_type>
struct box_interior_aug {
	using box_type = typename base_type::box_type;

	box_interior_aug() : box(base_type::get_empty_box())
	{
	}
	box_interior_aug(box_type const &b) : box(b)
	{
	}

	box_type const &get_box() const
	{
		return box;
	}

	template <typename leaf_type, typename interior_type>
	static box_type create(node *l, node *r)
	{
		return base_type::get_box(
			base_type::template retrieve_box<leaf_type,
							 interior_type>(l),
			base_type::template retrieve_box<leaf_type,
							 interior_type>(r));
	}

	template <typename leaf_type, typename interior_type>
	void update(node *l, node *r)
	{
		box = create<leaf_type, interior_type>(l, r);
	}

	template <typename leaf_type, typename interior_type,
		  typename TreeNodes>
	static box_type create(TreeNodes const &nodes)
	{
		box_type b = base_type::get_empty_box();
		for (auto *t : nodes)
			b = base_type::get_box(
				b, base_type::template retrieve_box<
					   leaf_type, interior_type>(t));
		return b;
	}

	template <typename leaf_type, typename interior_type,
		  typename TreeNodes>
	void update(TreeNodes const &nodes)
	{
		box = create<leaf_type, interior_type, TreeNodes>(nodes);
	}

	/*
	 * Unset means decide from the subtree size. Set overrides it, which is
	 * what batch delete needs: it shrinks sizes before the rebuild pass, so
	 * a size test at that point would under-parallelise.
	 */
	void set_parallel_flag(bool const flag)
	{
		force_par.emplace(flag);
	}

	void reset_parallel_flag()
	{
		force_par.reset();
	}

	bool get_parallel_flag_ini_status() const
	{
		return force_par.has_value();
	}

	bool force_parallel(size_t sz) const
	{
		return force_par.has_value()
			       ? force_par.value()
			       : sz > base_type::serial_build_cutoff;
	}

	void reset()
	{
		force_par.reset();
		box = base_type::get_empty_box();
	}

	box_type box;
	std::optional<bool> force_par;
};

} /* namespace psi */

#endif /* PSI_DEPENDENCE_AUGMENTATION_H_ */
