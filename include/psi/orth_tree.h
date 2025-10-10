#ifndef PSI_ORTH_TREE_H_
#define PSI_ORTH_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_based/orth_tree.h"

// Backward compatibility: alias pointer_based types into psi namespace
namespace psi {
using pointer_based::OrthTree;
}  // namespace psi

#endif  // PSI_ORTH_TREE_H_