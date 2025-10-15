#ifndef PSI_KD_TREE_H_
#define PSI_KD_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation

#include "array_view/kd_tree_array.h"
#include "pointer_view/kd_tree.h"

// Backward compatibility: alias pointer_view types into psi namespace
namespace psi {
using array_view::KdTreeArray;
using pointer_view::KdTree;
}  // namespace psi

#endif  // PSI_KD_TREE_H_
