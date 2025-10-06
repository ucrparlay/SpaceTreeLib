#ifndef PSI_P_TREE_IMPL_P_INTER_NODE_HPP_
#define PSI_P_TREE_IMPL_P_INTER_NODE_HPP_

#define PTREE_TEMPLATE template <typename Point, typename SplitRule, \
    uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
#define PTREE_CLASS PTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "../p_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi {}  // namespace psi

#undef PTREE_TEMPLATE
#undef PTREE_CLASS

#endif
