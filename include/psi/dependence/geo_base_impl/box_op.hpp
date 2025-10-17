#ifndef PSI_DEPENDENCE_GEO_BASE_IMPL_BOX_OP_HPP_
#define PSI_DEPENDENCE_GEO_BASE_IMPL_BOX_OP_HPP_

#include "dependence/geo_base.h"

namespace psi {

template <class TypeTrait>
inline typename GeoBase<TypeTrait>::Coord GeoBase<TypeTrait>::GetBoxMid(
    DimsType const d, Box const& box) {
  // TODO: change to Num::divide2ceil
  if constexpr (std::is_integral_v<Coord>) {
    // NOTE: since the points on the box line will be put in right, therefore
    // the box should always be rounded up as well.
    // Consideing the example (1, 1) and (1, 2), in order to split the points,
    // the new box should have split at y=2
    return static_cast<Coord>(
        std::ceil(static_cast<double>(box.first.pnt[d] + box.second.pnt[d]) /
                  static_cast<double>(2)));
  } else {
    return (box.first.pnt[d] + box.second.pnt[d]) / 2;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::LegalBox(Box const& box) {
  // TODO: remove it
  if (box == GetEmptyBox()) return true;

  if constexpr (kDim == 2) {
    return !Num::Gt(box.first.pnt[0], box.second.pnt[0]) &&
           !Num::Gt(box.first.pnt[1], box.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Gt(box.first.pnt[0], box.second.pnt[0]) &&
           !Num::Gt(box.first.pnt[1], box.second.pnt[1]) &&
           !Num::Gt(box.first.pnt[2], box.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Gt(box.first.pnt[i], box.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::WithinBox(Box const& a, Box const& b) {
  constexpr uint_fast8_t kDim = Point::GetDim();
  assert(LegalBox(a));
  assert(LegalBox(b));

  if constexpr (kDim == 2) {
    return !Num::Lt(a.first.pnt[0], b.first.pnt[0]) &&
           !Num::Lt(a.first.pnt[1], b.first.pnt[1]) &&
           !Num::Gt(a.second.pnt[0], b.second.pnt[0]) &&
           !Num::Gt(a.second.pnt[1], b.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Lt(a.first.pnt[0], b.first.pnt[0]) &&
           !Num::Lt(a.first.pnt[1], b.first.pnt[1]) &&
           !Num::Lt(a.first.pnt[2], b.first.pnt[2]) &&
           !Num::Gt(a.second.pnt[0], b.second.pnt[0]) &&
           !Num::Gt(a.second.pnt[1], b.second.pnt[1]) &&
           !Num::Gt(a.second.pnt[2], b.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(a.first.pnt[i], b.first.pnt[i]) ||
          Num::Gt(a.second.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::SameBox(Box const& a, Box const& b) {
  assert(LegalBox(a));
  assert(LegalBox(b));

  if constexpr (kDim == 2) {
    return Num::Eq(a.first.pnt[0], b.first.pnt[0]) &&
           Num::Eq(a.first.pnt[1], b.first.pnt[1]) &&
           Num::Eq(a.second.pnt[0], b.second.pnt[0]) &&
           Num::Eq(a.second.pnt[1], b.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return Num::Eq(a.first.pnt[0], b.first.pnt[0]) &&
           Num::Eq(a.first.pnt[1], b.first.pnt[1]) &&
           Num::Eq(a.first.pnt[2], b.first.pnt[2]) &&
           Num::Eq(a.second.pnt[0], b.second.pnt[0]) &&
           Num::Eq(a.second.pnt[1], b.second.pnt[1]) &&
           Num::Eq(a.second.pnt[2], b.second.pnt[2]);
    for (DimsType i = 0; i < kDim; ++i) {
      if (!Num::Eq(a.first.pnt[i], b.first.pnt[i]) ||
          !Num::Eq(a.second.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::WithinBox(Point const& p, Box const& box) {
  assert(LegalBox(box));

  if constexpr (kDim == 2) {
    return !Num::Lt(p.pnt[0], box.first.pnt[0]) &&
           !Num::Lt(p.pnt[1], box.first.pnt[1]) &&
           !Num::Gt(p.pnt[0], box.second.pnt[0]) &&
           !Num::Gt(p.pnt[1], box.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Lt(p.pnt[0], box.first.pnt[0]) &&
           !Num::Lt(p.pnt[1], box.first.pnt[1]) &&
           !Num::Lt(p.pnt[2], box.first.pnt[2]) &&
           !Num::Gt(p.pnt[0], box.second.pnt[0]) &&
           !Num::Gt(p.pnt[1], box.second.pnt[1]) &&
           !Num::Gt(p.pnt[2], box.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(p.pnt[i], box.first.pnt[i]) ||
          Num::Gt(p.pnt[i], box.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::BoxIntersectBox(Box const& a, Box const& b) {
  if (!LegalBox(a)) {
    std::cout << a.first << " " << a.second << std::endl;
  }
  assert(LegalBox(a) && LegalBox(b));

  if constexpr (kDim == 2) {
    return !(Num::Lt(a.second.pnt[0], b.first.pnt[0]) ||
             Num::Gt(a.first.pnt[0], b.second.pnt[0]) ||
             Num::Lt(a.second.pnt[1], b.first.pnt[1]) ||
             Num::Gt(a.first.pnt[1], b.second.pnt[1]));
  } else if constexpr (kDim == 3) {
    return !(Num::Lt(a.second.pnt[0], b.first.pnt[0]) ||
             Num::Gt(a.first.pnt[0], b.second.pnt[0]) ||
             Num::Lt(a.second.pnt[1], b.first.pnt[1]) ||
             Num::Gt(a.first.pnt[1], b.second.pnt[1]) ||
             Num::Lt(a.second.pnt[2], b.first.pnt[2]) ||
             Num::Gt(a.first.pnt[2], b.second.pnt[2]));
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(a.second.pnt[i], b.first.pnt[i]) ||
          Num::Gt(a.first.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::IsBoxLineInDimension(Box const& box,
                                                     DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(box.first[d], box.second[d]);
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineSplitBox(Coord const& line,
                                                     Box const& box,
                                                     DimsType d) {
  assert(LegalBox(box));
  return VerticalLineIntersectBoxExclude(line, box, d) ||
         (VerticalLineOnBoxRightEdge(line, box, d) &&
          !IsBoxLineInDimension(box, d));
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineOnBoxLeftEdge(Coord const& line,
                                                          Box const& box,
                                                          DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(line, box.first.pnt[d]);
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineOnBoxRightEdge(Coord const& line,
                                                           Box const& box,
                                                           DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(line, box.second.pnt[d]);
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineOnBoxEdge(Coord const& line,
                                                      Box const& box,
                                                      DimsType d) {
  assert(LegalBox(box));
  return VerticalLineOnBoxLeftEdge(line, box, d) ||
         VerticalLineOnBoxRightEdge(line, box, d);
}

template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineIntersectBox(Coord const& line,
                                                         Box const& box,
                                                         DimsType d) {
  assert(LegalBox(box));
  return Num::Geq(line, box.first.pnt[d]) && Num::Leq(line, box.second.pnt[d]);
}

// NOTE: if the line @line is one the boundary of the box, then it will be
// considered as not intersect
template <class TypeTrait>
inline bool GeoBase<TypeTrait>::VerticalLineIntersectBoxExclude(
    Coord const& line, Box const& box, DimsType d) {
  assert(LegalBox(box));
  return Num::Gt(line, box.first.pnt[d]) && Num::Lt(line, box.second.pnt[d]);
}

template <class TypeTrait>
inline typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetEmptyBox() {
  return Box(Point(std::numeric_limits<Coord>::max()),
             Point(std::numeric_limits<Coord>::lowest()));
}

template <class TypeTrait>
typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetBox(Box const& x,
                                                            Box const& y) {
  return Box(x.first.MinCoords(y.first), x.second.MaxCoords(y.second));
}

template <class TypeTrait>
template <typename Range>
typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetBox(Range&& In)
  requires(parlay::is_random_access_range_v<Range>)
{
  if constexpr (std::same_as<
                    typename std::remove_reference_t<Range>::value_type, Box>) {
    return GetBoxFromBoxSeq(std::forward<Range>(In));
  } else {
    return GetBoxFromSlice(parlay::make_slice(std::forward<Range>(In)));
  }
}

template <class TypeTrait>
template <typename SliceType>
typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetBoxFromSlice(
    SliceType const V) {
  if (V.size() == 0) {
    return GetEmptyBox();
  } else {
    auto minmax = [&](Box const& x, Box const& y) {
      return Box(x.first.MinCoords(y.first), x.second.MaxCoords(y.second));
    };
    auto boxes = parlay::delayed_seq<Box>(V.size(), [&](size_t i) {
      return Box(V[i].GetCoords(), V[i].GetCoords());
    });
    return parlay::reduce(boxes, parlay::make_monoid(minmax, boxes[0]));
  }
}

// NOTE: this function omit the possibility that T contains the bounding box --
// it will always try to reduce a bounding box in the tree T
template <class TypeTrait>
template <typename Leaf, typename Interior>
typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetBox(Node* T) {
  if (T->size == 0) {
    return GetEmptyBox();
  }
  if (T->is_leaf) {
    Leaf* TL = static_cast<Leaf*>(T);
    if (TL->is_dummy) {
      return Box(TL->pts[0], TL->pts[0]);
    }
    return GetBox(TL->pts.cut(0, TL->size));
    // auto nb = GetBox(TL->pts.cut(0, TL->size));
    // assert(SameBox(nb, TL->GetBox()));
    // return nb;
  }
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (psi::pointer_view::IsBinaryNode<Interior>) {
    Box left_box, right_box;
    parlay::par_do_if(
        T->size > kSerialBuildCutoff,
        [&]() { left_box = GetBox<Leaf, Interior>(TI->left); },
        [&]() { right_box = GetBox<Leaf, Interior>(TI->right); });
    // Box const& left_box = GetBox<Leaf, Interior>(TI->left);
    // Box const& right_box = GetBox<Leaf, Interior>(TI->right);
    return GetBox(left_box, right_box);
  } else if constexpr (psi::pointer_view::IsMultiNode<Interior>) {
    BoxSeq return_box_seq(Interior::GetRegions());
    parlay::parallel_for(
        0, Interior::GetRegions(),
        [&](size_t i) {
          return_box_seq[i] = GetBox<Leaf, Interior>(TI->tree_nodes[i]);
        },
        TI->size > kSerialBuildCutoff ? 1 : Interior::GetRegions());
    // for (size_t i = 0; i < Interior::GetRegions(); i++) {
    //   return_box_seq[i] = GetBox<Leaf, Interior>(TI->tree_nodes[i]);
    // }
    return GetBox(return_box_seq);
  } else {
    assert(false);
  }
}

template <class TypeTrait>
typename GeoBase<TypeTrait>::Box GeoBase<TypeTrait>::GetBoxFromBoxSeq(
    BoxSeq const& box_seq) {
  return parlay::reduce(
      box_seq, parlay::make_monoid(
                   [&](Box const& x, Box const& y) { return GetBox(x, y); },
                   GetEmptyBox()));
  // Box box = GetEmptyBox();
  // for (auto const& b : box_seq) {
  //   box = GetBox(box, b);
  // }
  // return std::move(box);
}

template <class TypeTrait>
typename GeoBase<TypeTrait>::Point GeoBase<TypeTrait>::GetBoxCenter(
    Box const& box) {
  Point center;
  if constexpr (kDim == 2) {
    center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
    center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
  } else if constexpr (kDim == 3) {
    center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
    center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
    center.pnt[2] = (box.first.pnt[2] + box.second.pnt[2]) / 2;
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      center.pnt[i] = GetBoxMid(i, box);
    }
  }
  return center;
}
}  // namespace psi

#endif  // PSI_DEPENDENCE_GEO_BASE_IMPL_BOX_OP_HPP_
