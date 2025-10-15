#ifndef PSI_BASE_TREE_H_
#define PSI_BASE_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_view/base_tree.h"

// Backward compatibility: alias pointer_view types into psi namespace
namespace psi {
using pointer_view::BaseTree;
}  // namespace psi

#endif  // PSI_BASE_TREE_H_
