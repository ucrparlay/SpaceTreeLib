#pragma once

#include "comparator.h"
#include "basic_point.h"

#include "baseTree.h"

#include "tree_node.h"
#include "utility/box_op.hpp"
#include "utility/dimensinality.hpp"
#include "utility/validation.hpp"
#include "utility/random_support.hpp"

#include "batch_op/build_tree.hpp"
// #include "batch_op/batch_insert.hpp"
// #include "batch_op/batch_delete.hpp"
// #include "batch_op/inner_tree.hpp"
#include "batch_op/batch_helpers.hpp"

#include "query_op/nn_search_helpers.h"
#include "query_op/nn_search.hpp"
#include "query_op/range_count.hpp"
#include "query_op/range_query.hpp"

#include "oct_tree.h"
#include "octTree/oct_build_tree.hpp"
#include "octTree/oct_inter_node.hpp"
