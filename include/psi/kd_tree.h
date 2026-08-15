#ifndef PSI_KD_TREE_H
#define PSI_KD_TREE_H

#include <array>
#include <functional>
#include <optional>
#include <utility>

#include "base_tree.h"
#include "dependence/augmentation.h"
#include "dependence/concepts.h"

namespace psi
{

template <typename Point, typename SplitRule,
	  typename LeafAugType = BoxLeafAug<BaseTree<Point>>,
	  typename InteriorAugType = BoxInteriorAug<BaseTree<Point>>,
	  uint_fast8_t kSkHeight = 6, uint_fast8_t kImbaRatio = 30>
class KdTree : public BaseTree<Point,
			       KdTree<Point, SplitRule, LeafAugType,
				      InteriorAugType, kSkHeight, kImbaRatio>,
			       kSkHeight, kImbaRatio>
{
public:
	using BT = BaseTree<Point,
			    KdTree<Point, SplitRule, LeafAugType,
				   InteriorAugType, kSkHeight, kImbaRatio>,
			    kSkHeight, kImbaRatio>;

	static_assert(
		LeafAugmentation<LeafAugType, typename BT::Slice>,
		"LeafAugType needs A(), A(Slice), UpdateAug(Slice), Reset()");
	static_assert(InteriorAugmentation<InteriorAugType>,
		      "InteriorAugType needs SetParallelFlag(bool), "
		      "ResetParallelFlag(), GetParallelFlagIniStatus(), "
		      "ForceParallel(size_t)");

	using BucketType = typename BT::BucketType;
	using BallsType = typename BT::BallsType;
	using DimsType = typename BT::DimsType;
	using BucketSeq = typename BT::BucketSeq;
	using Coord = typename Point::Coord;
	using Coords = typename Point::Coords;
	using Num = Num_Comparator<Coord>;
	using Slice = typename BT::Slice;
	using Points = typename BT::Points;
	using PointsIter = typename BT::PointsIter;
	using Box = typename BT::Box;
	using BoxSeq = typename BT::BoxSeq;

	using HyperPlane = typename BT::HyperPlane;
	using HyperPlaneSeq = typename BT::HyperPlaneSeq;
	using NodeBox = typename BT::NodeBox;
	using NodeBoxSeq = typename BT::NodeBoxSeq;
	using Splitter = HyperPlane;
	using SplitterSeq = HyperPlaneSeq;
	using SplitRuleType = SplitRule;

	static constexpr DimsType const kMD = 2;

	using Leaf = LeafNode<Point, Slice, BT::kLeafCapacity, LeafAugType,
			      parlay::move_assign_tag>;
	struct KdInteriorNode;
	using Interior = KdInteriorNode;

	using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
	using BoxCut = typename BT::BoxCut;

	template <typename Leaf, typename Interior, bool granularity,
		  typename... Args>
	friend Node *BT::RebuildSingleTree(Node *T, Args &&...args);

	template <typename Leaf, typename Interior, typename PrepareFunc,
		  typename... Args>
	friend Node *BT::RebuildWithInsert(Node *T, PrepareFunc prepare_func,
					   Slice In, Args &&...args);

	/* The split rule calls BuildRecursive back; see DivideSpace. */
	friend SplitRule;

	void KdTreeTag();

	template <typename Range>
	void Build(Range &&In);

	template <typename Range>
	void BatchInsert(Range &&In);

	// NOTE: every point is assumed to be in the tree; if that may not
	// hold, use BatchDiff
	template <typename Range>
	void BatchDelete(Range &&In);

	// NOTE: tolerates points that are not in the tree
	template <typename Range>
	void BatchDiff(Range &&In);

	template <typename Range>
	void Flatten(Range &&Out);

	template <typename Range>
	auto KNN(Node *T, Point const &q, kBoundedQueue<Point, Range> &bq);

	auto RangeCount(Box const &query_box);

	template <typename Range>
	auto RangeQuery(Box const &query_box, Range &&Out);

	constexpr void DeleteTree() override;

	constexpr static char const *GetTreeName()
	{
		return "KdTree";
	}

	constexpr static char const *CheckHasBox()
	{
		if constexpr (HasBox<InteriorAugType>)
			return "HasBox";
		else
			return "NoBox";
	}

private:
	void Build_(Slice In);
	Node *BuildRecursive(Slice In, Slice Out, DimsType dim, Box const &bx);
	Node *SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
				   Box const &bx);
	void DivideRotate(Slice In, SplitterSeq &pivots, DimsType dim,
			  BucketType idx, BoxSeq &box_seq, Box const &bx);
	void PickPivots(Slice In, size_t const &n, SplitterSeq &pivots,
			DimsType const dim, BoxSeq &box_seq, Box const &bx);

	void BatchInsert_(Slice In);
	Node *BatchInsertRecursive(Node *T, Slice In, Slice Out, DimsType d);

	void BatchDelete_(Slice In);
	NodeBox BatchDeleteRecursive(Node *T, Box const &bx, Slice In,
				     Slice Out, DimsType d, bool has_tomb);

	void BatchDiff_(Slice In);
	NodeBox BatchDiffRecursive(Node *T, Box const &bx, Slice In, Slice Out,
				   DimsType d);

	SplitRule split_rule_;
};

} // namespace psi

#include "kd_tree_impl/kd_batch_delete.hpp"
#include "kd_tree_impl/kd_batch_diff.hpp"
#include "kd_tree_impl/kd_batch_insert.hpp"
#include "kd_tree_impl/kd_build_tree.hpp"
#include "kd_tree_impl/kd_inter_node.hpp"
#include "kd_tree_impl/kd_override.hpp"

#endif // PSI_KD_TREE_H
