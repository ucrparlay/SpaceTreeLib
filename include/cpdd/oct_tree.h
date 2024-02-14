#pragma once

#include "comparator.h"
#include "baseTree.h"
#include "basic_point.h"
#include "query_op/nn_search_helpers.h"

namespace cpdd {

#define LOG  std::cout
#define ENDL std::endl << std::flush

template<typename point>
class octTree : public baseTree<point> {};

}  // namespace cpdd
