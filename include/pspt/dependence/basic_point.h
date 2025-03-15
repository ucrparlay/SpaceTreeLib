#ifndef PSPT_DEPENDENCE_BASIC_POINT_H_
#define PSPT_DEPENDENCE_BASIC_POINT_H_

#include <cstdint>
#include <tuple>
#include <variant>

#include "comparator.h"
#include "cpam/cpam.h"
#include "parlay/alloc.h"
#include "parlay/delayed.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"

namespace pspt {

template <typename T, uint_fast8_t d>
struct BasicPoint {
  using Coord = T;
  using Coords = std::array<T, d>;
  using Num = Num_Comparator<Coord>;
  using DimsType = uint_fast8_t;

  BasicPoint() {}

  explicit BasicPoint(T const val) { this->pnt.fill(val); }

  explicit BasicPoint(Coords const& _pnt) : pnt(_pnt) {}

  explicit BasicPoint(parlay::slice<T*, T*> x) {
    assert(x.size() == d);
    for (DimsType i = 0; i < d; ++i) {
      this->pnt[i] = x[i];
    }
  }

  explicit BasicPoint(T const* x) {
    for (DimsType i = 0; i < d; ++i) {
      this->pnt[i] = x[i];
    }
  }

  BasicPoint const MinCoords(BasicPoint const& b) const {
    BasicPoint p;

    if constexpr (d == 2) {
      p.pnt[0] = Num::Min(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::Min(this->pnt[1], b.pnt[1]);
    } else if constexpr (d == 3) {
      p.pnt[0] = Num::Min(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::Min(this->pnt[1], b.pnt[1]);
      p.pnt[2] = Num::Min(this->pnt[2], b.pnt[2]);
    } else {
      for (DimsType i = 0; i < d; ++i) {
        p.pnt[i] = Num::Min(this->pnt[i], b.pnt[i]);
      }
    }

    return p;
  }

  BasicPoint const MaxCoords(BasicPoint const& b) const {
    BasicPoint p;

    if constexpr (d == 2) {
      p.pnt[0] = Num::Max(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::Max(this->pnt[1], b.pnt[1]);
    } else if constexpr (d == 3) {
      p.pnt[0] = Num::Max(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::Max(this->pnt[1], b.pnt[1]);
      p.pnt[2] = Num::Max(this->pnt[2], b.pnt[2]);
    } else {
      for (DimsType i = 0; i < d; ++i) {
        p.pnt[i] = Num::Max(this->pnt[i], b.pnt[i]);
      }
    }

    return p;
  }

  static consteval auto GetDim() { return d; }

  bool SameDimension(BasicPoint const& b) const { return *this == b; }

  bool operator==(BasicPoint const& x) const {
    for (DimsType i = 0; i < d; ++i) {
      if (!Num::Eq(this->pnt[i], x.pnt[i])) return false;
    }
    return true;
  }

  bool operator<(BasicPoint const& x) const {
    for (DimsType i = 0; i < d; ++i) {
      if (Num::Lt(this->pnt[i], x.pnt[i]))
        return true;
      else if (Num::Gt(this->pnt[i], x.pnt[i]))
        return false;
      else
        continue;
    }
    return false;
  }

  Coord& operator[](DimsType i) { return pnt[i]; }

  Coord const& operator[](DimsType i) const { return pnt[i]; }

  friend std::ostream& operator<<(std::ostream& o, BasicPoint const& a) {
    o << "(";
    for (DimsType i = 0; i < d; ++i) {
      o << a.pnt[i] << (i == d - 1 ? "" : ", ");
    }
    o << ") " << std::flush;
    return o;
  }

  Coords pnt;
};

template <typename T, uint_fast8_t d, class AugType = std::monostate>
  requires std::totally_ordered<AugType>
struct AugPoint : BasicPoint<T, d> {
  using BT = BasicPoint<T, d>;
  using DimsType = BT::DimsType;

  AugPoint() {}

  explicit AugPoint(T const val) : BasicPoint<T, d>(val) {}

  explicit AugPoint(typename BasicPoint<T, d>::Coords const& _pnt)
      : BasicPoint<T, d>(_pnt) {}

  explicit AugPoint(parlay::slice<T*, T*> x) : BasicPoint<T, d>(x) {}

  explicit AugPoint(T const* x) : BasicPoint<T, d>(x) {}

  AugPoint(typename BasicPoint<T, d>::Coords const& _pnt, AugType const& aug)
      : BasicPoint<T, d>(_pnt), aug(aug) {}

  AugPoint(parlay::slice<T*, T*> x, AugType const& aug)
      : BasicPoint<T, d>(x), aug(aug) {}

  AugPoint(T const* x, AugType const& aug) : BasicPoint<T, d>(x), aug(aug) {}

  AugPoint(AugPoint const& p) : BasicPoint<T, d>(p), aug(p.aug) {}

  AugPoint(AugPoint&& p) noexcept
      : BasicPoint<T, d>(std::move(p)), aug(std::move(p.aug)) {}

  bool operator==(AugPoint const& x) const {
    return BasicPoint<T, d>::operator==(x) && aug == x.aug;
  }

  bool operator<(AugPoint const& x) const {
    return BasicPoint<T, d>::operator<(x)
               ? true
               : (BasicPoint<T, d>::operator==(x) ? aug < x.aug : false);
  }

  friend std::ostream& operator<<(std::ostream& o, AugPoint const& a) {
    o << "[" << BT::operator<<(o, a) << "-" << a.aug << "] " << std::flush;
    return o;
  }

  AugType aug;
};

}  // namespace pspt

#endif  // PSPT_DEPENDENCE_BASIC_POINT_H_
