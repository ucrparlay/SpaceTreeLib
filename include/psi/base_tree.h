#ifndef PSI_BASE_TREE_H_
#define PSI_BASE_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_based/base_tree.h"

// Backward compatibility: alias pointer_based types into psi namespace
namespace psi {
using pointer_based::BaseTree;
}  // namespace psi

#endif  // PSI_BASE_TREE_H_