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

  PointType() {}

  explicit PointType(T const val) { this->pnt.fill(val); }

  explicit PointType(Coords const& _pnt) : pnt(_pnt) {}

  explicit PointType(parlay::slice<T*, T*> x) {
    assert(x.size() == d);
    for (int i = 0; i < d; i++) {
      this->pnt[i] = x[i];
    }
  }

  explicit PointType(T const* x) {
    for (int i = 0; i < d; i++) {
      this->pnt[i] = x[i];
    }
  }

  inline PointType const MinCoords(PointType const& b) const {
    PointType p;
    for (uint_fast8_t i = 0; i < d; i++) {
      p.pnt[i] = Num::min(this->pnt[i], b.pnt[i]);
    }
    return std::move(p);
  }

  inline PointType const MaxCoords(PointType const& b) const {
    PointType p;
    for (uint_fast8_t i = 0; i < d; i++) {
      p.pnt[i] = Num::max(this->pnt[i], b.pnt[i]);
    }
    return std::move(p);
  }

  constexpr auto GetDim() const { return std::tuple_size_v<Coords>; }

  inline bool SameDimension(PointType const& b) const { return *this == b; }

  inline bool operator==(PointType const& x) const {
    for (uint_fast8_t i = 0; i < d; i++) {
      if (!Num::Eq(this->pnt[i], x.pnt[i])) return false;
    }
    return true;
  }

  inline bool operator<(PointType const& x) const {
    for (int i = 0; i < d; i++) {
      if (Num::Lt(this->pnt[i], x.pnt[i]))
        return true;
      else if (Num::Gt(this->pnt[i], x.pnt[i]))
        return false;
      else
        continue;
    }
    return false;
  }

  friend std::ostream& operator<<(std::ostream& o, PointType const& a) {
    o << "(";
    for (int i = 0; i < d; i++) {
      o << a.pnt[i] << (i == d - 1 ? "" : ", ");
    }
    o << ") " << std::flush;
    return o;
  }

  Coords pnt;
};

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

  inline PointID const MinCoords(PointID const& b) const {
    Coords pts;
    for (int i = 0; i < d; i++) {
      pts[i] = Num::min(this->pnt[i], b.pnt[i]);
    }
    return std::move(PointID(pts));
  }

  inline PointID const MaxCoords(PointID const& b) const {
    Coords pts;
    for (int i = 0; i < d; i++) {
      pts[i] = Num::max(this->pnt[i], b.pnt[i]);
    }
    return std::move(PointID(pts));
  }

  inline bool operator==(PointID const& x) const {
    for (int i = 0; i < d; i++) {
      if (!Num::Eq(this->pnt[i], x.pnt[i])) return false;
    }
    return this->id == x.id;
  }

  inline bool operator<(PointID const& x) const {
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
