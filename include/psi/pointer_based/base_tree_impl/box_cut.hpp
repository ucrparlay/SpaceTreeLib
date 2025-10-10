#ifndef PSI_POINTER_BASED_BASE_TREE_IMPL_BOX_CUT_HPP_
#define PSI_POINTER_BASED_BASE_TREE_IMPL_BOX_CUT_HPP_

#include "../base_tree.h"

namespace psi {
namespace pointer_based {
template <class TypeTrait, typename DerivedTree>
struct BaseTree<TypeTrait, DerivedTree>::BoxCut {
  using BT = BaseTree<TypeTrait, DerivedTree>;

  BoxCut(Box const& box, HyperPlane const& hp, bool go_left)
      : box(box), hp(hp), go_left(go_left) {}

  inline Box const& GetFirstBoxCut() {
    mod_dim = go_left ? &box.second.pnt[hp.second] : &box.first.pnt[hp.second];
    std::ranges::swap(hp.first, *mod_dim);
    return box;
  }

  inline Box const& GetSecondBoxCut() {
    std::ranges::swap(hp.first, *mod_dim);
    mod_dim = go_left ? &box.first.pnt[hp.second] : &box.second.pnt[hp.second];
    *mod_dim = hp.first;
    return box;
  }

  inline Box const& GetBox() const { return box; }

  inline HyperPlane const& GetHyperPlane() const { return hp; }

  Box box;
  Coord* mod_dim;
  HyperPlane hp;  // PARA: the split and the cutting dimension
  bool const go_left;
};
}  // namespace pointer_based
}  // namespace psi

#endif  // PSI_POINTER_BASED_BASE_TREE_IMPL_BOX_CUT_HPP_