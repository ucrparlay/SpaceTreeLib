#ifndef PSI_DEPENDENCE_AUGMENTATION_H_
#define PSI_DEPENDENCE_AUGMENTATION_H_

#include <cstddef>
#include <optional>

#include "tree_node.h"

namespace psi
{

/*
 * Default augmentations. KdTree and OrthTree use these unless told otherwise,
 * so the common case needs no user code.
 *
 * BT is a BaseTree used only for its box vocabulary, so BaseTree<Point> works
 * even though the tree instantiates BaseTree<Point, TheTree, ...>. Do not
 * reach for anything here that varies with the other parameters: BucketType
 * changes at kSkHeight > 7.
 */

/* Bounding box of the points in a leaf. */
template <typename BT>
struct BoxLeafAug {
	using Box = typename BT::Box;
	using Slice = typename BT::Slice;

	BoxLeafAug() : box(BT::GetEmptyBox())
	{
	}
	BoxLeafAug(Slice In) : box(BT::GetBox(In))
	{
	}

	Box const &GetBox() const
	{
		return box;
	}

	void UpdateAug(Slice In)
	{
		box = BT::GetBox(In);
	}

	void Reset()
	{
		box = BT::GetEmptyBox();
	}

	Box box;
};

/*
 * Bounding box of a subtree, plus the flag that picks parallel granularity.
 *
 * Create and Update are member templates because they are handed the very
 * Interior type that holds this augmentation. Both a binary and a multi-way
 * form are present; the node calls Create<Leaf, Interior>(...) and the last
 * parameter is deduced from the argument.
 */
template <typename BT>
struct BoxInteriorAug {
	using Box = typename BT::Box;

	BoxInteriorAug() : box(BT::GetEmptyBox())
	{
	}
	BoxInteriorAug(Box const &b) : box(b)
	{
	}

	Box const &GetBox() const
	{
		return box;
	}

	template <typename Leaf, typename Interior>
	static Box Create(Node *l, Node *r)
	{
		return BT::GetBox(BT::template RetrieveBox<Leaf, Interior>(l),
				  BT::template RetrieveBox<Leaf, Interior>(r));
	}

	template <typename Leaf, typename Interior>
	void Update(Node *l, Node *r)
	{
		box = Create<Leaf, Interior>(l, r);
	}

	template <typename Leaf, typename Interior, typename TreeNodes>
	static Box Create(TreeNodes const &nodes)
	{
		Box b = BT::GetEmptyBox();
		for (auto *t : nodes)
			b = BT::GetBox(
				b, BT::template RetrieveBox<Leaf, Interior>(t));
		return b;
	}

	template <typename Leaf, typename Interior, typename TreeNodes>
	void Update(TreeNodes const &nodes)
	{
		box = Create<Leaf, Interior, TreeNodes>(nodes);
	}

	/*
	 * Unset means decide from the subtree size. Set overrides it, which is
	 * what batch delete needs: it shrinks sizes before the rebuild pass, so
	 * a size test at that point would under-parallelise.
	 */
	void SetParallelFlag(bool const flag)
	{
		force_par.emplace(flag);
	}

	void ResetParallelFlag()
	{
		force_par.reset();
	}

	bool GetParallelFlagIniStatus() const
	{
		return force_par.has_value();
	}

	bool ForceParallel(size_t sz) const
	{
		return force_par.has_value() ? force_par.value()
					     : sz > BT::kSerialBuildCutoff;
	}

	void Reset()
	{
		force_par.reset();
		box = BT::GetEmptyBox();
	}

	Box box;
	std::optional<bool> force_par;
};

} // namespace psi

#endif // PSI_DEPENDENCE_AUGMENTATION_H_
