#pragma once
#include <parlay/internal/get_time.h>

#include "cpam_sample_sort.h"
#include "devotail_integer_sort.h"
#include "get_time.h"
#include "parlay/internal/integer_sort.h"
#include "parlay/primitives.h"

namespace cpam {

template <class Entry>
struct build {
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using ET = typename Entry::entry_t;
  using filling_curve_t = typename Entry::filling_curve_t;
  using sort_output_value_t = typename Entry::sort_output_value_t;

  constexpr static auto less = [](const ET& a, const ET& b) {
    return Entry::comp(Entry::get_key(a), Entry::get_key(b));
  };

  // sorts a sequence, then removes all but first element with equal keys
  // the sort is not necessarily stable, so any element could be kept
  template <class Seq>
  // static parlay::sequence<ET> sort_remove_duplicates(
  static auto sort_remove_duplicates(Seq const& A) {  // ?? const
    // BUG: should add the handling of the empty sequence
    // if (A.size() == 0) return parlay::sequence<ET>(0);

    // parlay::internal::timer t("");
    // assert(parlay::all_of(
    //     A, [&](auto const& p) { return std::get<0>(p).first == 0; }));
    // auto B = parlay::internal::cpam::cpam_sample_sort<filling_curve_t>(
    //     parlay::make_slice(A.begin(), A.end()),
    //     [](auto const& a, auto const& b) {
    //       return std::get<0>(a).first < std::get<0>(b).first;
    //     });
    auto B = parlay::internal::cpam::cpam_sample_sort<filling_curve_t,
                                                      sort_output_value_t>(
        parlay::make_slice(A.begin(), A.end()),
        [&](auto const& a, auto const& b) { return a.first < b.first; });
    // auto B = parlay::internal::cpam::cpam_sample_sort<filling_curve_t>(
    //     parlay::make_slice(A.begin(), A.end()), less);
    // auto B = parlay::internal::integer_sort(
    //     parlay::make_slice(A.begin(), A.end()),
    //     [](auto const& k) { return Entry::get_key(k).code; });
    // auto B = parlay::cpam::integer_sort2(
    //     parlay::make_slice(A.begin(), A.end()),
    //     [](auto const& k) { return Entry::get_key(k).code; });
    // t.next("sort");

    // auto Fl = parlay::delayed_seq<bool>(
    //     B.size(), [&](size_t i) { return (i == 0) || less(B[i - 1], B[i]);
    //     });
    //
    // auto o = parlay::pack(B, Fl);
    // t.next("pack");
    // return o;
    return B;
  }

  template <class Seq, class Reduce>
  static parlay::sequence<ET> sort_reduce_duplicates(Seq const& A,
                                                     Reduce const& reduce) {
    using E = typename Seq::value_type;
    using Vi = typename E::second_type;

    timer t("sort_reduce_duplicates", false);
    size_t n = A.size();
    if (n == 0) return parlay::sequence<ET>(0);
    auto lessE = [](E const& a, E const& b) {
      return Entry::comp(a.first, b.first);
    };

    auto B = parlay::internal::sample_sort(parlay::make_slice(A), lessE);
    t.next("sort");

    // determines the index of start of each block of equal keys
    // and copies values into vals
    parlay::sequence<bool> Fl(n);
    auto Vals = parlay::tabulate(n, [&](size_t i) -> Vi {
      Fl[i] = (i == 0) || lessE(B[i - 1], B[i]);
      return B[i].second;
    });
    t.next("copy, set flags");

    auto I = parlay::pack_index<node_size_t>(Fl);
    t.next("pack index");

    // combines over each block of equal keys using function reduce
    auto a = parlay::tabulate(I.size(), [&](size_t i) -> ET {
      size_t start = I[i];
      size_t end = (i == I.size() - 1) ? n : I[i + 1];
      return ET(B[start].first, reduce(Vals.cut(start, end)));
    });
    t.next("reductions");
    // tabulate set over all entries of i
    return a;
  }

  template <class Seq, class Reduce>
  static parlay::sequence<ET> sort_reduce_duplicates_a(Seq const& A,
                                                       Reduce const& reduce) {
    using E = typename Seq::value_type;
    using Vi = typename E::second_type;

    timer t("sort_reduce_duplicates", false);
    size_t n = A.size();
    if (n == 0) return parlay::sequence<ET>(0);
    auto lessE = [](E const& a, E const& b) {
      return Entry::comp(a.first, b.first);
    };

    auto B = parlay::internal::sample_sort(parlay::make_slice(A), lessE);
    t.next("sort");

    // determines the index of start of each block of equal keys
    // and copies values into vals
    parlay::sequence<bool> Fl(n);  // ?? should it be unitialized
    auto Vals = parlay::tabulate(n, [&](size_t i) -> Vi {
      Fl[i] = (i == 0) || lessE(B[i - 1], B[i]);
      return B[i].second;
    });
    t.next("copy, set flags");

    auto I = parlay::pack_index<node_size_t>(Fl);
    t.next("pack index");

    // combines over each block of equal keys using function reduce
    auto a = parlay::tabulate(I.size(), [&](size_t i) -> ET {
      size_t start = I[i];
      size_t end = (i == I.size() - 1) ? n : I[i + 1];
      return ET(B[start].first, reduce(Vals.cut(start, end)));
    });
    t.next("reductions");
    // tabulate set over all entries of i
    return a;
  }

  template <class Seq, class Bin_Op>
  static parlay::sequence<ET> sort_combine_duplicates(Seq const& A, Bin_Op& f) {
    auto mon = parlay::make_monoid(f, V());
    auto reduce_op = [&](parlay::slice<V*, V*> S) {  // ?? should it be &
      return parlay::reduce(S, mon);
    };
    return sort_reduce_duplicates_a(A, reduce_op);
  }

  template <class Seq, class Bin_Op>
  static parlay::slice<ET*, ET*> sort_combine_duplicates_inplace(Seq const& A,
                                                                 Bin_Op& f) {
    auto less = [&](ET a, ET b) { return Entry::comp(a.first, b.first); };
    parlay::internal::quicksort(A.begin(), A.size(), less);
    size_t j = 0;
    for (size_t i = 1; i < A.size(); i++) {
      if (less(A[j], A[i]))
        A[++j] = A[i];
      else
        A[j].second = f(A[j].second, A[i].second);
    }
    return A.cut(0, j + 1);
  }
};

}  // namespace cpam
