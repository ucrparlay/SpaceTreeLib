#ifndef PSTP_DEPENDENCE_SEARCH_CONTAINER_H_
#define PSTP_DEPENDENCE_SEARCH_CONTAINER_H_

#include "comparator.h"
#include "parlay/slice.h"

namespace pstp {

template <typename Point, typename T>
class NN_Comparator {
  using Coord = typename Point::Coord;
  using Num = Num_Comparator<Coord>;

 public:
  bool operator()(T const& a, T const& b) {
    return Num::Lt(a.second, b.second);
  }
};

template <typename Point, typename T,
          typename Compare = NN_Comparator<Point, T>>
class kBoundedQueue {
  // NOTE: A priority queue with fixed maximum capacity. While the queue has
  // not reached its maximum capacity, elements are inserted as they will be
  // in a heap, the root (top()) being such that Compare(top(),x)=false for
  // any x in the queue. Once the queue is full, trying to insert x in the
  // queue will have no effect if Compare(x,top())=false. Otherwise, the
  // element at the root of the heap is removed and x is inserted so as to
  // keep the heap property.

  // NOTE:
  // A simplified version of CGAL bounded_priority_queue
  // https://github.com/CGAL/cgal/blob/v5.4/Spatial_searching/include/CGAL/Spatial_searching/internal/bounded_priority_queue.h

  using Coord = typename Point::Coord;
  // using T = std::pair<Point*, Coord>;

 public:
  kBoundedQueue(Compare const& comp = Compare()) : m_comp(comp) {}

  kBoundedQueue(parlay::slice<T*, T*> data_slice,
                Compare const& comp = Compare())
      : m_count(0), m_data(data_slice), m_comp(comp) {}

  /** Sets the max number of elements in the queue */
  void resize(parlay::slice<T*, T*> data_slice) {
    m_data = data_slice;
    m_count = 0;
  }

  void reset() { m_count = 0; }

  inline bool full() const { return m_count == m_data.size(); }

  inline T const& top() const { return m_data[0]; }

  inline Coord const top_value() const { return m_data[0].second; }

  inline parlay::slice<T*, T*> data() const { return m_data; }

  inline void insert(std::pair<std::reference_wrapper<Point>, Coord> const x) {
    // T x( _x.first, _x.second );
    T* data1 = (&m_data[0] - 1);
    if (full()) {
      if (m_comp(x, top())) {
        // insert x in the heap at the correct place,
        // going down in the tree.
        size_t j(1), k(2);
        while (k <= m_count) {
          T* z = &(data1[k]);
          if ((k < m_count) && m_comp(*z, data1[k + 1])) z = &(data1[++k]);
          if (m_comp(*z, x)) break;
          data1[j] = *z;
          j = k;
          k = j << 1;  // a son of j in the tree
        }
        data1[j] = x;
      }
    } else {
      // insert element as in a heap
      size_t i(++m_count), j(0);
      while (i >= 2) {
        j = i >> 1;  // father of i in the tree
        T& y = data1[j];
        if (m_comp(x, y)) break;
        data1[i] = y;
        i = j;
      }
      data1[i] = x;
    }
  }

 public:
  size_t m_count = 0;
  parlay::slice<T*, T*> m_data;
  Compare m_comp;
};

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_SEARCH_CONTAINER_H_
