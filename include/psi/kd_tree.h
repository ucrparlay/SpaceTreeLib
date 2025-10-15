#ifndef PSI_KD_TREE_H_
#define PSI_KD_TREE_H_

// Backward compatibility header
// This file forwards to the pointer-based implementation

#include "array_based/kd_tree_array.h"
#include "pointer_based/kd_tree.h"

// Backward compatibility: alias pointer_based types into psi namespace
namespace psi {
using array_based::KdTreeArray;
using pointer_based::KdTree;

}  // namespace psi

#endif  // PSI_KD_TREE_H_
