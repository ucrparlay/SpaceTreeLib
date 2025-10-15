#ifndef PSI_COVER_TREE_H_
#define PSI_COVER_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation
#include "pointer_view/cover_tree.h"

// Backward compatibility: alias pointer_view types into psi namespace
namespace psi {
using pointer_view::CoverTree;
}  // namespace psi

#endif  // PSI_COVER_TREE_H_
