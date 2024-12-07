#pragma once

#include <cstdint>
#include <tuple>

#include "comparator.h"
#include "parlay/alloc.h"
#include "parlay/delayed.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"

namespace cpdd {

template <typename T, uint_fast8_t d>
struct PointType {
  using Coord = T;
  using Coords = std::array<T, d>;
  using Num = Num_Comparator<Coord>;
  using DimsType = uint_fast8_t;

  PointType() {}

  explicit PointType(T const val) { this->pnt.fill(val); }

  explicit PointType(Coords const& _pnt) : pnt(_pnt) {}

  explicit PointType(parlay::slice<T*, T*> x) {
    assert(x.size() == d);
    for (DimsType i = 0; i < d; i++) {
      this->pnt[i] = x[i];
    }
  }

  explicit PointType(T const* x) {
    for (DimsType i = 0; i < d; i++) {
      this->pnt[i] = x[i];
    }
  }

  PointType const MinCoords(PointType const& b) const {
    PointType p;

    if constexpr (d == 2) {
      p.pnt[0] = Num::min(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::min(this->pnt[1], b.pnt[1]);
    } else if constexpr (d == 3) {
      p.pnt[0] = Num::min(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::min(this->pnt[1], b.pnt[1]);
      p.pnt[2] = Num::min(this->pnt[2], b.pnt[2]);
    } else {
      for (DimsType i = 0; i < d; i++) {
        p.pnt[i] = Num::min(this->pnt[i], b.pnt[i]);
      }
    }

    return p;
  }

  PointType const MaxCoords(PointType const& b) const {
    PointType p;

    if constexpr (d == 2) {
      p.pnt[0] = Num::max(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::max(this->pnt[1], b.pnt[1]);
    } else if constexpr (d == 3) {
      p.pnt[0] = Num::max(this->pnt[0], b.pnt[0]);
      p.pnt[1] = Num::max(this->pnt[1], b.pnt[1]);
      p.pnt[2] = Num::max(this->pnt[2], b.pnt[2]);
    } else {
      for (DimsType i = 0; i < d; i++) {
        p.pnt[i] = Num::max(this->pnt[i], b.pnt[i]);
      }
    }

    return p;
  }

  static consteval auto GetDim() { return d; }

  bool SameDimension(PointType const& b) const { return *this == b; }

  bool operator==(PointType const& x) const {
    for (DimsType i = 0; i < d; i++) {
      if (!Num::Eq(this->pnt[i], x.pnt[i])) return false;
    }
    return true;
  }

  bool operator<(PointType const& x) const {
    for (DimsType i = 0; i < d; i++) {
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

  friend std::ostream& operator<<(std::ostream& o, PointType const& a) {
    o << "(";
    for (DimsType i = 0; i < d; i++) {
      o << a.pnt[i] << (i == d - 1 ? "" : ", ");
    }
    o << ") " << std::flush;
    return o;
  }

  Coords pnt;
};

// NOTE: sample usage of point with id associated
template <typename T, uint_fast8_t d, typename IDtype = uint_fast64_t>
struct PointID : PointType<T, d> {
  using Coord = T;
  using Coords = std::array<T, d>;
  using Num = Num_Comparator<Coord>;
  using ID = IDtype;

  PointID() {}

  explicit PointID(T const val) : id(0) { this->pnt.fill(val); }

  explicit PointID(Coords const& _pnt) : PointType<T, d>(_pnt), id(0) {}

  PointID(Coords const& _pnt, ID _id) : PointType<T, d>(_pnt), id(_id) {}

  PointID(parlay::slice<T*, T*> x, ID _id) : PointType<T, d>(x), id(_id) {}

  PointID(T const* x, ID _id) : PointType<T, d>(x), id(_id) {}

  PointID const MinCoords(PointID const& b) const {
    Coords pts;
    for (int i = 0; i < d; i++) {
      pts[i] = Num::min(this->pnt[i], b.pnt[i]);
    }
    return std::move(PointID(pts));
  }

  PointID const MaxCoords(PointID const& b) const {
    Coords pts;
    for (int i = 0; i < d; i++) {
      pts[i] = Num::max(this->pnt[i], b.pnt[i]);
    }
    return std::move(PointID(pts));
  }

  bool operator==(PointID const& x) const {
    for (int i = 0; i < d; i++) {
      if (!Num::Eq(this->pnt[i], x.pnt[i])) return false;
    }
    return this->id == x.id;
  }

  bool operator<(PointID const& x) const {
    if (this->id == x.id) {
      for (int i = 0; i < d; i++) {
        if (Num::Lt(this->pnt[i], x.pnt[i]))
          return true;
        else if (Num::Gt(this->pnt[i], x.pnt[i]))
          return false;
        else
          continue;
      }
      return false;
    } else {
      return this->id < x.id;
    }
  }

  friend std::ostream& operator<<(std::ostream& o, PointID const& a) {
    o << a.id << "-";
    o << "(";
    for (int i = 0; i < d; i++) {
      o << a.pnt[i] << ", ";
    }
    o << ") " << std::flush;
    return o;
  }

  ID get_id() { return id; }

  ID id;
};

}  // namespace cpdd
