#ifndef PSI_ORTH_TREE_H_
#define PSI_ORTH_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_view/orth_tree.h"

// Backward compatibility: alias pointer_view types into psi namespace
namespace psi {
using pointer_view::OrthTree;
}  // namespace psi

#endif  // PSI_ORTH_TREE_H_
