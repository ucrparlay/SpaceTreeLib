/*
 * A p_tree past 2^32 points counts them all. Its CPAM nodes are made by hand:
 * two children stand in for 2^31 points each.
 */

#include <cstddef>
#include <utility>

#include "psi/p_tree.h"
#include "tests/unit/check.h"

int main()
{
	using node_type =
		psi::cpam::basic_node<psi::cpam::weight_balanced_tree::data,
				      std::pair<long, long>,
				      psi::cpam::null_encoder, 40>;
	CASE("p/size-past-2^32");

	auto *left = node_type::single({1, 1});
	auto *right = node_type::single({3, 3});
	auto *root = node_type::single({2, 2});
	left->s = size_t(1) << 31;
	right->s = size_t(1) << 31;
	root->lc = left;
	root->rc = right;
	node_type::update(root);
	CHECK_EQ(size_t(node_type::size(root)), (size_t(1) << 32) + 1);

	node_type::free_node(left);
	node_type::free_node(right);
	node_type::free_node(root);
	return psi_test_result();
}
