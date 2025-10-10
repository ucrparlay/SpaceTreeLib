#ifndef PSI_P_TREE_H_
#define PSI_P_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_based/p_tree.h"

// Backward compatibility: alias pointer_based types into psi namespace
namespace psi {
using pointer_based::PTree;
}  // namespace psi

#endif  // PSI_P_TREE_H_