#pragma once
#include "basic_node_helpers.h"
#include "parlay/internal/binary_search.h"
#include "parlay/primitives.h"
#include "utils.h"

// *******************************************
//   MAPS and SETS
// *******************************************

namespace cpam {

template <class Seq, class EntryT>
struct map_ops : Seq {
  using Entry = EntryT;
  using node = typename Seq::node;
  using regular_node = typename Seq::regular_node;
  using compressed_node = typename Seq::compressed_node;
  using ET = typename Seq::ET;
  using GC = typename Seq::GC;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using _Seq = Seq;
  using ptr = typename GC::ptr;
  using Seq::B;
  using Seq::kBaseCaseSize;

  static bool comp(K a, K b) { return Entry::comp(a, b); }
  // static bool comp(K const& a, K const& b) { return Entry::comp(a, b); }
  static K get_key(node* s) { return Entry::get_key(*Seq::get_entry_p(s)); }
  static V get_val(node* s) { return Entry::get_val(*Seq::get_entry_p(s)); }
  // static K get_key(node *s) { return Entry::get_key(Seq::get_entry(s));}
  // static V get_val(node *s) { return Entry::get_val(Seq::get_entry(s));}

  static std::optional<ET> find_compressed(ptr b, K const& key) {
    std::optional<ET> ret;
    auto find_f = [&](const ET& et) {
      if (!(Entry::comp(key, Entry::get_key(et)) ||
            Entry::comp(Entry::get_key(et), key))) {
        ret = et;
      }
    };
    Seq::iterate_seq(b.node_ptr(), find_f);
    // TODO: should probably use decode_cond or something like this.
    return ret;
  }

  static std::optional<ET> find(ptr b, K const& key) {
    if (b.empty()) return {};
    if (b.is_compressed()) {
      return find_compressed(std::move(b), key);
    }
    auto [lc, e, rc, m] = Seq::expose(std::move(b));
    std::optional<ET> ret = {};
    if (Entry::comp(key, Entry::get_key(e))) {
      ret = find(std::move(lc), key);
    } else if (Entry::comp(Entry::get_key(e), key)) {
      ret = find(std::move(rc), key);
    } else {
      ret = {e};
    }
    GC::decrement(m);
    return ret;
  }

  static std::optional<ET> find_compressed2(node* b, K const& key) {
    std::optional<ET> ret;
    auto find_f = [&](const ET& et) -> bool {
      if (!(Entry::comp(key, Entry::get_key(et)) ||
            Entry::comp(Entry::get_key(et), key))) {
        ret = et;
        return false;
      }
      return true;
    };
    Seq::iterate_cond(b, find_f);
    return ret;
  }

  static std::optional<ET> find2(node* b, K const& key) {
    if (!b) return {};
    // if (Seq::is_compressed(b)) return find_compressed2(b, key);
    auto f = [&](const ET& et) { return Entry::get_key(et); };
    if (Seq::is_compressed(b))
      return Seq::find_compressed(b, f, Entry::comp, key);

    auto [lc, e, rc] = Seq::expose_simple(b);
    if (Entry::comp(key, Entry::get_key(e)))
      return find2(lc, key);
    else if (Entry::comp(Entry::get_key(e), key))
      return find2(rc, key);
    else
      return {e};
  }

  template <class BinaryOp>
  static inline void update_value(regular_node* a, BinaryOp const& op) {
    ET re = Seq::get_entry(a);
    Entry::set_val(re, op(re));
    Seq::set_entry(a, re);
  }

  // TODO
  //  static node* previous(node* b, const K& key) {
  //    node* r = NULL;
  //    while (b) {
  //      if (Entry::comp(get_key(b), key)) {r = b; b = b->rc;}
  //      else b = b->lc;
  //    }
  //    return r;
  //  }
  //
  //  static node* next(node* b, const K& key) {
  //    node* r = NULL;
  //    while (b) {
  //      if (Entry::comp(key, get_key(b)) ) {r = b; b = b->lc;}
  //      else b = b->rc;
  //    }
  //    return r;
  //  }

  static size_t rank(node* b, K const& key, size_t sum = 0) {
    if (!b) return sum;
    if (Seq::is_compressed(b)) {
      auto fn = [&](auto const& et) {
        if (Entry::comp(Entry::get_key(et), key)) {
          sum++;
        }
      };
      Seq::iterate_seq(b, fn);
      return sum;
    }
    auto rb = Seq::cast_to_regular(b);
    if (Entry::comp(get_key(rb), key)) {
      return rank(rb->rc, key, sum + Seq::size(rb->lc) + 1);
    } else if (!Entry::comp(key, get_key(rb))) {  // equal
      return sum + Seq::size(rb->lc);             // +1?
    }
    return rank(rb->lc, key, sum);
  }

  struct split_info {
    node* l;
    std::optional<ET> mid;
    node* r;
    split_info(node* l, std::optional<ET> mid, node* r)
        : l(l), mid(mid), r(r) {};
  };

  static split_info split(ptr a, K const& k) {
    if (a.empty()) return split_info(NULL, std::nullopt, NULL);
    auto [lc, e, rc, m] = Seq::expose(std::move(a));
    K const& kmid = Entry::get_key(e);

    if (Entry::comp(kmid, k)) {
      split_info bstpair = split(std::move(rc), k);
      bstpair.l = Seq::join(lc.node_ptr(), e, bstpair.l, m);
      return bstpair;
    } else if (Entry::comp(k, kmid)) {
      split_info bstpair = split(std::move(lc), k);
      bstpair.r = Seq::join(bstpair.r, e, rc.node_ptr(), m);
      return bstpair;
    } else {
      GC::decrement(m);
      return split_info(lc.node_ptr(), e, rc.node_ptr());
    }
  }

  // TODO
  //  static split_info split_copy(node* bst, const K& e) {
  //    if (!bst) return split_info(NULL, NULL, false);
  //
  //    else if (Entry::comp(get_key(bst), e)) {
  //      // following should be a copy ????
  //      node* join = Seq::make_node(Seq::get_entry(bst));
  //      split_info bstpair = split_copy(bst->rc, e);
  //      GC::increment(bst->lc);
  //      bstpair.first = Seq::node_join(bst->lc, bstpair.first, join);
  //      return bstpair;
  //    }
  //    else if (Entry::comp(e, get_key(bst))) {
  //      node* join = Seq::make_node(Seq::get_entry(bst));
  //      split_info bstpair = split_copy(bst->lc, e);
  //      GC::increment(bst->rc);
  //      bstpair.second = Seq::node_join(bstpair.second, bst->rc, join);
  //      return bstpair;
  //    }
  //    else {
  //      GC::increment(bst->lc); GC::increment(bst->rc);
  //      split_info ret(bst->lc, bst->rc, true);
  //      ret.entry = Seq::get_entry(bst);
  //      return ret;
  //    }
  //  }

  // TODO: deprecate?
  //  // A version that will reuse a node if ref count is 1
  //  // Will decrement ref count of root if it is copied
  //  static split_info split_inplace(node* bst, const K& e) {
  //    if (!bst) return split_info(NULL, NULL, false);
  //    else if (bst->ref_cnt > 1) {
  //		//std::cout << "rc>1" << std::endl;
  //      split_info ret = split_copy(bst, e);
  //      GC::decrement_recursive(bst);
  //      return ret;
  //    }
  //    else if (Entry::comp(get_key(bst), e)) {
  //      split_info bstpair = split_inplace(bst->rc, e);
  //      bstpair.first = Seq::node_join(bst->lc, bstpair.first, bst);
  //      return bstpair;
  //    }
  //    else if (Entry::comp(e, get_key(bst))) {
  //      split_info bstpair = split_inplace(bst->lc, e);
  //      bstpair.second = Seq::node_join(bstpair.second, bst->rc, bst);
  //      return bstpair;
  //    }
  //    else {
  //      split_info ret(bst->lc, bst->rc, true);
  //      ret.entry = Seq::get_entry(bst);
  //      GC::decrement(bst);
  //      return ret;
  //    }
  //  }
  //
  //  static inline split_info split(node* bst, const K& e) {
  //    return split_inplace(bst, e);
  //  }

  template <class BinaryOp>
  static inline void combine_values(ET& re, ET e, bool reverse,
                                    BinaryOp const& op) {
    if (reverse)
      Entry::set_val(re, op(Entry::get_val(e), Entry::get_val(re)));
    else
      Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(e)));
  }

  // a is the root, e is the new (incoming value), rev=false
  // combine(old, new)
  template <class BinaryOp>
  static inline void combine_values(node* a, ET e, bool reverse,
                                    BinaryOp const& op) {
    //    ET& re = Seq::get_entry(a);
    //    combine_values(re, e, reverse, op);
    ET re = Seq::get_entry(a);
    if (reverse)
      Entry::set_val(re, op(Entry::get_val(e), Entry::get_val(re)));
    else
      Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(e)));
    Seq::set_entry(a, re);
  }

  // Value-based version of combine. No reverse on this one for some
  // reason.
  template <class VE, class BinaryOp>
  static inline void combine_values_v(node* a, VE v0, BinaryOp const& op) {
    ET re = Seq::get_entry(a);
    Entry::set_val(re, op(Entry::get_val(re), v0));
    Seq::set_entry(a, re);
  }

  template <class VE, class BinaryOp>
  static inline void update_valuev(node* a, VE v0, BinaryOp const& op) {
    ET re = Seq::get_entry(a);
    Entry::set_val(re, op(Entry::get_val(re), v0));
    Seq::set_entry(a, re);
  }

  // valid only if key type is printable
  static void output_r(node* t) {
    if (!t) return;
    output_r(t->lc);
    std::cout << get_key(t) << " ";
    output_r(t->rc);
  }

  // valid only if key type is printable
  static void output(node* t) {
    output_r(t);
    std::cout << std::endl;
  }

  // reuses node
  template <class BinaryOp>
  static node* uniont(ptr b1, ptr b2, BinaryOp const& op) {
    if (b1.empty()) return b2.node_ptr();
    if (b2.empty()) return b1.node_ptr();

    size_t n1 = b1.size(), n2 = b2.size();
    if (n1 + n2 <= kBaseCaseSize) {
      return union_bc(std::move(b1), std::move(b2), op);
    }

    auto [l2, e, r2, m2] = Seq::expose(std::move(b2));
    if (!m2) m2 = Seq::single(e);

    auto sp1 = split(std::move(b1), Entry::get_key(e));
#ifdef DEBUG
    Seq::check_structure(sp1.l);
    Seq::check_structure(sp1.r);
#endif

    auto [l, r] = utils::fork<node*>(
        Seq::do_parallel(n1, n2),
        [&]() { return uniont(sp1.l, std::move(l2), op); },
        [&]() { return uniont(sp1.r, std::move(r2), op); });
#ifdef DEBUG
    Seq::check_structure(l);
    Seq::check_structure(r);
#endif

    if (sp1.mid) combine_values(m2, *sp1.mid, true, op);
    return Seq::node_join(l, r, m2);
  }

  template <class Seq1, class Seq2, class BinaryOp>
  static node* intersect(typename Seq1::node* b1, typename Seq2::ptr b2,
                         BinaryOp const& op) {
    if (!b1) return NULL;
    if (b2.empty()) {
      Seq1::GC::decrement_recursive(b1);
      return NULL;
    }

    size_t n1 = Seq1::size(b1);
    size_t n2 = b2.size();

    if (n1 + n2 <= kBaseCaseSize) {
      return intersect_bc<Seq1, Seq2, BinaryOp>(std::move(b1), std::move(b2),
                                                op);
    }

    auto [l2, e2, r2, m2] = Seq2::expose(std::move(b2));
    auto key = Seq2::Entry::get_key(e2);
    auto sp1 = Seq1::split(b1, key);

#ifdef DEBUG
    Seq::check_structure(sp1.l);
    Seq::check_structure(sp1.r);
#endif

    auto [l, r] = utils::fork<node*>(
        Seq::do_parallel(n1, n2),
        [&]() { return intersect<Seq1, Seq2>(sp1.l, std::move(l2), op); },
        [&]() { return intersect<Seq1, Seq2>(sp1.r, std::move(r2), op); });

#ifdef DEBUG
    Seq::check_structure(l);
    Seq::check_structure(r);
#endif

    Seq2::GC::decrement(m2);
    if (sp1.mid) {
      ET e(key, op(Seq1::Entry::get_val(*sp1.mid), Seq2::Entry::get_val(e2)));
      return Seq::join(l, e, r, nullptr);
    } else
      return Seq::join2(l, r);
  }

  static node* difference(ptr b1, node* b2) {
    if (b1.empty()) {
      GC::decrement_recursive(b2);
      return NULL;
    }
    if (!b2) return b1.node_ptr();

    size_t n1 = b1.size(), n2 = Seq::size(b2);
    if (n1 + n2 <= kBaseCaseSize) {
      return difference_bc(std::move(b1), std::move(b2));
    }

    auto [l1, e, r1, m1] = Seq::expose(std::move(b1));
    auto sp2 = split(b2, Entry::get_key(e));

    auto [l, r] = utils::fork<node*>(
        Seq::do_parallel(n1, n2),
        [&]() { return difference(std::move(l1), sp2.l); },
        [&]() { return difference(std::move(r1), sp2.r); });

    if (sp2.mid) {
      GC::decrement(m1);
      return Seq::join2(l, r);
    } else {
      if (!m1) m1 = Seq::single(e);
      return Seq::node_join(l, r, m1);
    }
  }

  // TODO: this seems inefficient.
  static node* range_root(ptr b, K const& key_left, K const& key_right) {
    if (b.empty()) return NULL;
    if (b.is_compressed()) {
      return b.node_ptr();
    }
    auto k = Entry::get_key(b.entry());
    if (Entry::comp(key_right, k)) {
      auto [lc, e, rc, m] = Seq::expose(std::move(b));
      auto ret = range_root(std::move(lc), key_left, key_right);
      GC::decrement(m);
      return ret;
    } else if (Entry::comp(k, key_left)) {
      auto [lc, e, rc, m] = Seq::expose(std::move(b));
      auto ret = range_root(std::move(rc), key_left, key_right);
      GC::decrement(m);
      return ret;
    }
    return b.node_ptr();
  }

  static node* range_root_2(node* b, K const& key_left, K const& key_right) {
    if (!b) return NULL;
    if (Seq::is_compressed(b)) {
      return b;
    }
    auto k = Entry::get_key(Seq::get_entry(b));
    auto rb = Seq::cast_to_regular(b);
    if (Entry::comp(key_right, k)) {
      return range_root_2(rb->lc, key_left, key_right);
    } else if (Entry::comp(k, key_left)) {
      return range_root_2(rb->rc, key_left, key_right);
    }
    return b;
  }

  template <class F>
  static node* comp_bc(ptr b, F& f) {
    auto r = b.node_ptr();

    ET stack[2 * B];
    size_t offset = 0;
    auto copy_f = [&](ET a) {
      auto k = Entry::get_key(a);
      if (f(k)) {
        parlay::move_uninitialized(stack[offset++], a);
      }
    };
    Seq::iterate_seq(r, copy_f);
    assert(offset <= 2 * B);

    if (offset < B) {
      return Seq::to_tree_impl((ET*)stack, offset);
    } else {
      return Seq::make_compressed(stack, offset);
    }
  }

  static node* left(ptr b, K const& e) {
    if (b.empty()) return NULL;
    if (b.is_compressed()) {
      auto comp = [&](K const& k) { return k <= e; };
      return comp_bc(std::move(b), comp);
    }
    auto [lc, b_e, rc, m] = Seq::expose(std::move(b));
    auto k = Entry::get_key(b_e);
    if (Entry::comp(e, k)) {
      GC::decrement(m);
      return left(std::move(lc), e);
    }
    return Seq::join(lc.node_ptr(), b_e, left(std::move(rc), e), m);
  }

  static node* left_2(node* b, K const& e) {
    if (!b) return NULL;
    if (Seq::is_compressed(b)) {
      auto comp = [&](K const& k) { return k <= e; };
      return comp_bc(b, comp);
    }
    auto [lc, b_e, rc] = Seq::expose_simple(b);
    auto k = Entry::get_key(b_e);
    if (Entry::comp(e, k)) {
      return left_2(lc, e);
    }
    GC::increment(lc);
    return Seq::join(lc, b_e, left_2(rc, e), nullptr);
  }

  static node* right(ptr b, K const& e) {
    if (b.empty()) return NULL;
    if (b.is_compressed()) {
      auto comp = [&](K const& k) { return e <= k; };
      return comp_bc(std::move(b), comp);
    }
    auto [lc, b_e, rc, m] = Seq::expose(std::move(b));
    auto k = Entry::get_key(b_e);
    if (Entry::comp(k, e)) {
      GC::decrement(m);
      return right(std::move(rc), e);
    }
    return Seq::join(right(std::move(lc), e), b_e, rc.node_ptr(), m);
  }

  // TODO: get rid of right() and replace with right_2(). Ditto with
  // left() and range() (and their _2 counterparts).
  static node* right_2(node* b, K const& e) {
    if (!b) return NULL;
    if (Seq::is_compressed(b)) {
      auto comp = [&](K const& k) { return e <= k; };
      return comp_bc(b, comp);
    }
    auto [lc, b_e, rc] = Seq::expose_simple(b);
    auto k = Entry::get_key(b_e);
    if (Entry::comp(k, e)) {
      return right_2(rc, e);
    }
    GC::increment(rc);
    return Seq::join(right_2(lc, e), b_e, rc, nullptr);
  }

  static node* range(ptr b, K const& low, K const& high) {
    node* r = range_root(std::move(b), low, high);
    if (!r) return NULL;
    if (Seq::is_compressed(r)) {
      // Build tree only on elements in the range in the leaf.
      auto comp = [&](K const& k) { return low <= k && k <= high; };
      return comp_bc(std::move(b), comp);
    } else {
      regular_node* rr = Seq::cast_to_regular(r);
      regular_node* root = Seq::single(Seq::get_entry(r));
      auto lc = rr->lc;
      auto rc = rr->rc;
      auto ret = Seq::node_join(right(ptr(lc, true), low),
                                left(ptr(rc, true), high), root);

      GC::decrement_recursive(r);
      return ret;
    }
  }

  static node* range_2(node* b, K const& low, K const& high) {
    node* r = range_root_2(b, low, high);
    if (!r) return NULL;
    if (Seq::is_compressed(r)) {
      // Build tree only on elements in the range in the leaf.
      auto comp = [&](K const& k) { return low <= k && k <= high; };
      return comp_bc(r, comp);
    } else {
      auto [lc, e, rc] = Seq::expose_simple(r);
      return Seq::join(right_2(lc, low), e, left_2(rc, high), nullptr);
    }
  }

  // TODO: move both of the following fns to a better location.
  // templatized version of inc_if
  template <bool copy = false>
  static node* inc_tmpl(node* x) {
    if constexpr (copy) {
      GC::increment(x);
    }
    return x;
  }
  template <bool copy = false>
  static regular_node* make_node_tmpl(regular_node* x) {
    if constexpr (copy) {
      return Seq::make_regular_node(Seq::get_entry(x));
    }
    return x;
  }

  template <class Func, class J>
  static node* insert_compressed(compressed_node* b, const ET& e,
                                 Func const& f) {
    ET arr[2 * B + 1];
    Seq::compressed_node_elms(b, arr);
    size_t n = Seq::size(b);
    assert(n <= 2 * B);
    ET merged[2 * B + 1];
    size_t out_off = 0;
    size_t k = 0;
    K key = Entry::get_key(e);
    bool placed = false;
    while (k < n) {
      if (Entry::comp(Entry::get_key(arr[k]), key)) {
        parlay::move_uninitialized(merged[out_off++], arr[k++]);
      } else if (Entry::comp(key, Entry::get_key(arr[k]))) {
        parlay::assign_uninitialized(merged[out_off++], e);
        placed = true;
        break;
      } else {  // arr[k] == key
        parlay::assign_uninitialized(merged[out_off], arr[k++]);
        combine_values(merged[out_off], e, true, f);
        out_off++;
        placed = true;
        break;
      }
    }
    while (k < n) {
      parlay::move_uninitialized(merged[out_off++], arr[k++]);
    }
    if (!placed) {
      parlay::assign_uninitialized(merged[out_off++], e);
    }
    return Seq::make_compressed(merged, out_off);
  }

  template <class Func>
  static node* insert_compressed2(node* b, const ET& e, Func const& f) {
    ET merged[2 * B + 1];

    K key = Entry::get_key(e);
    size_t out_off = 0;
    bool placed = false;
    auto merge = [&](const ET& et) {
      if (!placed) {
        if (Entry::comp(Entry::get_key(et), key)) {
          parlay::assign_uninitialized(merged[out_off++], et);
        } else if (Entry::comp(key, Entry::get_key(et))) {
          parlay::assign_uninitialized(merged[out_off++], e);
          placed = true;
        } else {  // get_key(et) == key
          parlay::assign_uninitialized(merged[out_off], et);
          combine_values(merged[out_off], e, true, f);
          out_off++;
          placed = true;
        }
      } else {
        parlay::assign_uninitialized(merged[out_off++], et);
      }
    };
    Seq::iterate_seq(b, merge);
    if (!placed) {
      parlay::assign_uninitialized(merged[out_off++], e);
    }
    return Seq::make_compressed(merged, out_off);
  }

  // A specialized version of insert that will switch to the version with
  // copy=true once it hits a node with ref_cnt(node) > 1.
  template <class Func, class J, bool copy = false>
  static node* insert_tmpl(node* b, const ET& e, Func const& f, J const& join) {
    if (!b) return Seq::single(e);

    if constexpr (!copy) {
      if (Seq::ref_cnt(b) > 1) {
        auto r = insert_tmpl<Func, J, true>(b, e, f, join);
        GC::decrement_recursive(b);
        return r;
      }
    }
    if (Seq::is_compressed(b)) {
      auto r = insert_compressed2<Func>(b, e, f);
      if constexpr (!copy) {
        GC::decrement_recursive(b);
      }
      return r;
    }

    regular_node* br = Seq::cast_to_regular(b);
    regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(br));
    if (Entry::comp(get_key(br), Entry::get_key(e))) {
      node* r = insert_tmpl<Func, J, copy>(br->rc, e, f, join);
      return join(inc_tmpl<copy>(br->lc), r, o);
    } else if (Entry::comp(Entry::get_key(e), get_key(br))) {
      node* l = insert_tmpl<Func, J, copy>(br->lc, e, f, join);
      return join(l, inc_tmpl<copy>(br->rc), o);
    } else {
      ET be = Seq::get_entry(br);
      Seq::set_entry(o, e);
      combine_values(o, be, true, f);
      return join(inc_tmpl<copy>(br->lc), inc_tmpl<copy>(br->rc), o);
    }
  }

  template <class Func>
  static node* insert(node* b, const ET& e, Func const& f) {
    auto join = [](node* l, node* r, regular_node* m) {
      auto ret = Seq::node_join(l, r, m);
      return ret;
    };
    return insert_tmpl<Func, decltype(join), false>(b, e, f, join);
  }

  //  static node* remove_compressed(node* nb, const K& key) {
  //    assert(Seq::is_compressed(nb));
  //    auto b = Seq::cast_to_compressed(nb);
  //    ET arr[2*B];
  //    Seq::compressed_node_elms(b, arr);
  //    size_t n = Seq::size(b);
  //    assert(n <= 2*B);
  //    ET merged[2*B];
  //    size_t k = 0;
  //    size_t out_off = 0;
  //    while (k < n) {
  //      if (Entry::comp(Entry::get_key(arr[k]), key)) {
  //        parlay::move_uninitialized(merged[out_off++], arr[k++]);
  //      } else if (Entry::comp(key, Entry::get_key(arr[k]))) {
  //        break;
  //      } else {  // arr[k] == key
  //        arr[k].~ET();
  //        k++;
  //        out_off++;
  //        break;
  //      }
  //    }
  //    while (k < n) {
  //      parlay::move_uninitialized(merged[out_off++], arr[k++]);}
  //    assert(out_off >= B);  // Is this always guaranteed? TODO: check
  //    return Seq::make_compressed(merged, out_off);
  //  }

  static node* remove_compressed2(node* nb, K const& key) {
    auto b = Seq::cast_to_compressed(nb);
    ET merged[2 * B];

    size_t out_off = 0;
    bool found = false;
    auto merge = [&](const ET& et) {
      if (!found) {
        if (Entry::comp(Entry::get_key(et), key)) {
          parlay::assign_uninitialized(merged[out_off++], et);
        } else if (Entry::comp(key, Entry::get_key(et))) {
          parlay::assign_uninitialized(merged[out_off++], et);
        } else {  // get_key(et) == key
          found = true;
        }
      } else {
        parlay::assign_uninitialized(merged[out_off++], et);
      }
    };
    Seq::iterate_seq(b, merge);

    if (out_off < B) {
      return Seq::to_tree_impl((ET*)merged, out_off);
    }
    return Seq::make_compressed(merged, out_off);
  }

  // A specialized version of insert that will switch to the version with
  // copy=true once it hits a node with ref_cnt(node) > 1.
  template <class J, bool copy = false>
  static node* remove_tmpl(node* b, K const& k, J const& join) {
    if (!b) return Seq::empty();

    if constexpr (!copy) {
      if (Seq::ref_cnt(b) > 1) {
        std::cout << "in remove: ref_cnt = " << Seq::ref_cnt(b) << std::endl;
        auto r = remove_tmpl<J, true>(b, k, join);
        GC::decrement_recursive(b);
        return r;
      }
    }

    if (Seq::is_compressed(b)) {
      auto r = remove_compressed2(b, k);
      if constexpr (!copy) {
        GC::decrement_recursive(b);
      }
      return r;
    }

    regular_node* br = Seq::cast_to_regular(b);
    if (Entry::comp(get_key(br), k)) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(br));
      node* r = remove_tmpl<J, copy>(br->rc, k, join);
      return join(inc_tmpl<copy>(br->lc), r, o);
    } else if (Entry::comp(k, get_key(br))) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(br));
      node* l = remove_tmpl<J, copy>(br->lc, k, join);
      return join(l, inc_tmpl<copy>(br->rc), o);
    } else {
      auto lc = br->lc, rc = br->rc;
      auto ret = Seq::join2(inc_tmpl<copy>(lc), inc_tmpl<copy>(rc));
      GC::decrement_recursive(br);
      std::cout << "Called join2!" << std::endl;
      return ret;
    }
  }

  //  static node* deletet(node* b, const K& k) {
  //    auto join = [] (node* l, node* r, regular_node* m) {
  //      auto ret = Seq::node_join(l,r,m);
  //      return ret;
  //    };
  //    return remove_tmpl<decltype(join), false>(b, k, join);
  //  }

  static node* remove_compressed3(ptr b, K const& k) {
    auto r = b.node_ptr();
    ET merged[2 * B];
    size_t out_off = 0;
    auto merge = [&](const ET& et) {
      if (Entry::comp(Entry::get_key(et), k) ||
          Entry::comp(k, Entry::get_key(et))) {
        parlay::assign_uninitialized(merged[out_off++], et);
      }
    };
    Seq::iterate_seq(r, merge);

    if (out_off < B) {
      return Seq::to_tree_impl((ET*)merged, out_off);
    }
    return Seq::make_compressed(merged, out_off);
  }

  template <class J>
  static node* my_remove_2(ptr b, K const& k, J const& join) {
    if (b.empty()) return nullptr;
    if (b.is_compressed()) {
      return remove_compressed3(std::move(b), k);
    }
    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    K ek = Entry::get_key(e);
    if (Entry::comp(ek, k)) {
      regular_node* o = (root != nullptr) ? root : Seq::single(e);
      auto ret = my_remove(std::move(rc), k, join);
      return join(lc.node_ptr(), ret, o);
    } else if (Entry::comp(k, ek)) {
      regular_node* o = (root != nullptr) ? root : Seq::single(e);
      auto ret = my_remove(std::move(lc), k, join);
      return join(ret, rc.node_ptr(), o);
    }
    GC::decrement(root);
    return Seq::join2(lc.node_ptr(), rc.node_ptr());
  }

  template <class J>
  static node* my_remove(ptr b, K const& k, J const& join) {
    if (b.empty()) return nullptr;
    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    K ek = Entry::get_key(e);
    if (Entry::comp(ek, k)) {
      regular_node* o = (root != nullptr) ? root : Seq::single(e);
      auto ret = my_remove(std::move(rc), k, join);
      return join(lc.node_ptr(), ret, o);
    } else if (Entry::comp(k, ek)) {
      regular_node* o = (root != nullptr) ? root : Seq::single(e);
      auto ret = my_remove(std::move(lc), k, join);
      return join(ret, rc.node_ptr(), o);
    }
    GC::decrement(root);
    return Seq::join2(lc.node_ptr(), rc.node_ptr());
  }

  static node* deletet(node* b, K const& k) {
    auto join = [](node* l, node* r, regular_node* m) {
      auto ret = Seq::node_join(l, r, m);
      return ret;
    };
    //     return my_remove<decltype(join)>(b, k, join);
    return my_remove_2<decltype(join)>(b, k, join);
    //    return remove_tmpl<decltype(join), false>(b, k, join);
  }

  template <class BinaryOp>
  static node* update_compressed(node* b, K const& key, BinaryOp const& op) {
    ET arr[2 * B];
    Seq::compressed_node_elms(b, arr);
    size_t n = Seq::size(b);
    assert(n <= 2 * B);
    ET merged[2 * B];
    size_t k = 0;
    size_t out_off = 0;
    while (k < n) {
      if (Entry::comp(Entry::get_key(arr[k]), key)) {
        parlay::move_uninitialized(merged[out_off++], arr[k++]);
      } else if (Entry::comp(key, Entry::get_key(arr[k]))) {
        break;
      } else {
        ET re = arr[k];
        Entry::set_val(re, op(re));
        parlay::move_uninitialized(merged[out_off++], re);
        k++;
        break;
      }
    }
    while (k < n) {
      parlay::move_uninitialized(merged[out_off++], arr[k++]);
    }
    return Seq::make_compressed(merged, out_off);
  }

  template <class Func, class J, bool copy = false>
  static node* update_tmpl(node* b, K const& k, Func const& f, J const& join) {
    if (!b) return b;
    if constexpr (!copy) {
      if (Seq::ref_cnt(b) > 1) {
        auto r = update_tmpl<Func, J, true>(b, k, f, join);
        GC::decrement_recursive(b);
        return r;
      }
    }

    if (Seq::is_compressed(b)) {
      auto r = update_compressed(b, k, f);
      if constexpr (!copy) {
        GC::decrement_recursive(b);
      }
      return r;
    }

    regular_node* br = Seq::cast_to_regular(b);
    if (Entry::comp(get_key(br), k)) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(br));
      node* r = update_tmpl<Func, J, copy>(br->rc, k, f, join);
      return join(inc_tmpl<copy>(br->lc), r, o);
    } else if (Entry::comp(k, get_key(b))) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(br));
      node* l = update_tmpl<Func, J, copy>(br->lc, k, f, join);
      return join(l, inc_tmpl<copy>(br->rc), o);
    } else {
      node* l = inc_tmpl<copy>(br->lc);
      node* r = inc_tmpl<copy>(br->rc);
      regular_node* o = make_node_tmpl<copy>(br);

      update_value(o, f);
      return join(l, r, o);
    }
  }

  template <class Func>
  static node* update(node* b, K const& k, Func const& f,
                      bool extra_ptr = false) {
    auto join = [](node* l, node* r, regular_node* m) {
      return Seq::node_join(l, r, m);
    };
    return update_tmpl<Func, decltype(join), false>(b, k, f, join);
  }

  static node* multidelete_bc(ptr b1, K* A, size_t n) {
    auto n_b1 = b1.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) { parlay::move_uninitialized(stack[offset++], a); };
    Seq::iterate_seq(n_b1, copy_f);

    // parlay::sort_inplace(
    //     parlay::make_slice(stack, stack + offset),
    //     [&](const ET& a, const ET& b) { return comp(a.aug, b.aug); });
    // parlay::sort_inplace(parlay::make_slice(A, A + n),
    //                      [&](K const& a, K const& b) { return comp(a, b); });

    Seq::decrement_recursive(n_b1);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = A[j];
      // std::cout << stack[i] << " " << k_b << std::endl;
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        // puts("xxxxxxxxxxxxxxxxxxxxxxxxx");
        j++;
      } else {  // equals, delete element.
        assert(k_a.id == k_b.id);
        assert(k_a.code == k_b.code);
        // std::cout << "delete " << k_a << "=" << k_b << std::endl;
        // puts("-------------------------");
        stack[i].~ET();
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  static node* multi_delete_sorted_size(ptr b, K* A, size_t n, size_t& sz) {
    if (n == 0) return b.node_ptr();
    sz++;
    if (b.empty()) return nullptr;

    // size_t tot = b.size() + n;
    if (b.size() <= kBaseCaseSize) {
      sz += count_size_ptr(std::move(b));
      sz--;
      return nullptr;
      // std::cout << "kBaseCaseSize" << std::endl;
      // return multidelete_bc(std::move(b), A, n);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));

    K bk = Entry::get_key(e);
    auto less_val = [&](K& a) -> bool { return Entry::comp(a, bk); };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));

    auto P = utils::fork<node*>(
        // true,  // Seq::do_parallel(b.size(), n),
        // !(mid == 0 || mid == n),
        false,  // Seq::do_parallel(b.size(), n),
        [&]() { return multi_delete_sorted_size(std::move(lc), A, mid, sz); },
        [&]() {
          return multi_delete_sorted_size(std::move(rc), A + mid + dup,
                                          n - mid - dup, sz);
        });

    if (!dup) {
      if (!root) root = Seq::single(e);
      return Seq::node_join(P.first, P.second, root);
    } else {
      GC::decrement(root);
      return Seq::join2(P.first, P.second);
    }
  }

  // assumes array A is of length n and is sorted with no duplicates
  static node* multi_delete_sorted(ptr b, K* A, size_t n) {
    if (b.empty()) return nullptr;
    if (n == 0) return b.node_ptr();

    // size_t tot = b.size() + n;
    if (b.size() <= kBaseCaseSize) {
      // std::cout << "kBaseCaseSize" << std::endl;
      return multidelete_bc(std::move(b), A, n);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));

    K bk = Entry::get_key(e);
    auto less_val = [&](K& a) -> bool { return Entry::comp(a, bk); };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));

    auto P = utils::fork<node*>(
        // true,  // Seq::do_parallel(b.size(), n),
        !(mid == 0 || mid == n),
        // false,  // Seq::do_parallel(b.size(), n),
        [&]() { return multi_delete_sorted(std::move(lc), A, mid); },
        [&]() {
          return multi_delete_sorted(std::move(rc), A + mid + dup,
                                     n - mid - dup);
        });

    if (!dup) {
      if (!root) root = Seq::single(e);
      return Seq::node_join(P.first, P.second, root);
    } else {
      GC::decrement(root);
      return Seq::join2(P.first, P.second);
    }
  }

  // TODO(laxman): This is more like a multi-update?
  template <class VE, class CombineOp>
  static node* multidelete_sorted_map_bc(ptr b1, std::pair<K, VE>* A, size_t n,
                                         CombineOp const& combine_op) {
    auto n_b1 = b1.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) { parlay::move_uninitialized(stack[offset++], a); };
    Seq::iterate_seq(n_b1, copy_f);
    assert(offset <= kBaseCaseSize);

    Seq::decrement_recursive(n_b1);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = A[j].first;
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, combine_op(Entry::get_val(re), A[j].second));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  // assumes array A is of length n and is sorted with no duplicates
  // combine: <V, VE> -> V
  // map: VE -> V
  template <class VE, class CombineOp>
  static node* multi_delete_sorted_map(ptr b, std::pair<K, VE>* A, size_t n,
                                       CombineOp const& combine_op) {
    if (b.empty()) return nullptr;
    if (n == 0) return b.node_ptr();

    // size_t tot = b.size() + n;
    // TODO: add a flag to control base-case granularity?
    if (b.size() <= kBaseCaseSize) {
      return multidelete_sorted_map_bc(std::move(b), A, n, combine_op);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    if (!root) root = Seq::single(e);

    K bk = Entry::get_key(e);
    auto less_val = [&](std::pair<K, VE>& a) -> bool {
      return Entry::comp(a.first, bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid].first));

    auto P = utils::fork<node*>(
        true,  // Seq::do_parallel(b.size(), n),
        [&]() {
          return multi_delete_sorted_map(std::move(lc), A, mid, combine_op);
        },
        [&]() {
          return multi_delete_sorted_map(std::move(rc), A + mid + dup,
                                         n - mid - dup, combine_op);
        });

    if (dup) combine_values_v(root, A[mid].second, combine_op);
    return Seq::node_join(P.first, P.second, root);
  }

  // assumes array A is of length n and is sorted with no duplicates
  // template <class BinaryOp>
  // static node* multi_insert_sorted(ptr b, ET* A, size_t n, BinaryOp const&
  // op) {
  // static size_t count_size(ptr b) {
  //   if (b.empty()) return 0;
  //   if (b.is_compressed()) {
  //     return 1;
  //   }
  //   auto [lc, e, rc, root] = Seq::expose(std::move(b));
  //   size_t l, r;
  //   parlay::par_do_if(
  //       root->s >= 100000, [&]() { l = count_size(std::move(lc)); },
  //       [&]() { r = count_size(std::move(rc)); });
  //   // l = count_size(std::move(lc));
  //   // r = count_size(std::move(rc));
  //   return l + r + 1;
  // }
  static size_t count_size(node* b) {
    if (!b) return 0;
    if (Seq::is_compressed(b)) {
      return 1;
    }
    // auto [lc, e, rc, root] = Seq::expose(std::move(b));
    // auto [lc, e, rc] = Seq::expose_simple(b);
    auto rb = Seq::cast_to_regular(b);
    size_t l, r;
    parlay::par_do_if(
        rb->s >= 1000, [&]() { l = count_size(rb->lc); },
        [&]() { r = count_size(rb->rc); });
    // l = count_size(rb->lc);
    // r = count_size(rb->rc);
    return l + r + 1;
  }

  static size_t count_size_ptr(ptr b) {
    if (b.empty()) return 0;
    if (b.is_compressed()) {
      return 1;
    }
    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    // auto [lc, e, rc] = Seq::expose_simple(b);
    // auto rb = Seq::cast_to_regular(b);
    size_t l, r;
    // parlay::par_do([&]() { l = count_size_ptr(std::move(lc)); },
    //                [&]() { r = count_size_ptr(std::move(rc)); });
    l = count_size_ptr(std::move(lc));
    r = count_size_ptr(std::move(rc));
    return l + r + 1;
  }

  template <class BinaryOp, typename T>
  static node* multiinsert_compress(ptr b1, T* A, size_t n,
                                    BinaryOp const& op) {
    auto n_b1 = b1.node_ptr();
    // ET stack[kBaseCaseSize + 1];
    // size_t offset = 0;
    // auto copy_f = [&](ET& a) {
    //   parlay::move_uninitialized(stack[offset++], a);
    // };
    // Seq::iterate_seq(n_b1, copy_f);
    // Seq::decrement_recursive(n_b1);
    size_t offset = Seq::size(n_b1);
    auto c = Seq::cast_to_compressed(n_b1);
    uint8_t* data_start =
        (((uint8_t*)c) + sizeof(typename Seq::aug_compressed_node));
    ET* stack = (ET*)data_start;

    ET output[offset + n + 1];
    // auto output = parlay::sequence<ET>::uninitialized(n + offset + 1);

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      // auto const& k_b = Entry::get_key(A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      auto const& k_b = Entry::get_key(et);
      if (comp(k_a, k_b)) {
        // parlay::move_uninitialized(output[out_off++], stack[i]);
        // parlay::assign_uninitialized(output[out_off++], stack[i]);
        output[out_off++] = stack[i];
        i++;
      }
      // else {
      //   // WARN: this assumes there is no duplicates
      //   // WARN: if is, the duplicates will be replaced automatically
      //   auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      //   output[out_off++] = et;
      //   j++;
      // }
      else if (comp(k_b, k_a)) {
        // parlay::move_uninitialized(output[out_off++], A[j]);
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        // parlay::move_uninitialized(output[out_off++], et);
        output[out_off++] = et;
        j++;
      } else {
        // parlay::move_uninitialized(output[out_off], stack[i]);
        // parlay::assign_uninitialized(output[out_off], stack[i]);
        output[out_off] = stack[i];
        ET& re = output[out_off];
        // Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(A[j])));
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(et)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      // parlay::move_uninitialized(output[out_off++], stack[i]);
      // parlay::assign_uninitialized(output[out_off++], stack[i]);
      output[out_off++] = stack[i];
      i++;
    }
    while (j < nB) {
      // parlay::move_uninitialized(output[out_off++], A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      // parlay::move_uninitialized(output[out_off++], et);
      output[out_off++] = et;
      j++;
    }

    // NOTE: the new returned tree replace the n_b1, decrease the counter
    Seq::decrement_recursive(n_b1);

    // build returned tree
    if (out_off < B) {
      // return Seq::to_tree_impl((ET*)output, out_off);
      return Seq::to_tree_impl(output, out_off);
    } else {
      // return Seq::make_compressed(output, out_off);
      return Seq::from_array(output, out_off);
    }
  }

  template <class BinaryOp, typename T>
  static node* multi_insert_sorted(ptr b, T* A, size_t n, BinaryOp const& op) {
    // std::cout << "multi_insert_sorted" << std::endl;
    if (n == 0) return b.node_ptr();

    if (b.empty()) {
      std::cout << "empty " << n << std::endl;
      node* x = Seq::from_array(A, n);
      return x;
    }

    size_t tot = b.size() + n;
    // if (tot <= kBaseCaseSize) {
    if (b.is_compressed()) {
      if (tot <= 2 * B) {  // can put in one leaf directly
        return multiinsert_bc_append(std::move(b), A, n, op);
      } else if (n <= 4 * B) {
        // return multiinsert_bc(std::move(b), A, n, op);
        return multiinsert_compress(std::move(b), A, n, op);
      } else {  // the merge overhead may be large, should expose compress
        ;
      }
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    if (!root) root = Seq::single(e);

    K bk = Entry::get_key(e);
    // auto less_val = [&](ET& a) -> bool {
    //   return Entry::comp(Entry::get_key(a), bk);
    // };
    auto less_val = [&](auto const& a) -> bool {
      auto et = basic_node_helpers::get_entry_indentity<ET>(a);
      return Entry::comp(Entry::get_key(et), bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    // size_t mid = n / 2;
    // std::cout << e << " " << mid << " " << n << std::endl;
    // bool dup = (mid < n) && (!Entry::comp(bk, Entry::get_key(A[mid])));
    bool dup =
        (mid < n) &&
        (!Entry::comp(
            bk, Entry::get_key(
                    basic_node_helpers::get_entry_indentity<ET>(A, mid))));

    auto P = utils::fork<node*>(
        // true,  // Seq::do_parallel(b.size(), n),
        n >= 1024 && !(mid == 0 || mid == n),
        // 0,
        // n >= 512,
        // Seq::do_parallel(b.size(), n),
        [&]() { return multi_insert_sorted(std::move(lc), A, mid, op); },
        [&]() {
          return multi_insert_sorted(std::move(rc), A + mid + dup,
                                     n - mid - dup, op);
        });

    // if (dup) combine_values(root, A[mid], false, op);
    if (dup) {
      auto mid_et = basic_node_helpers::get_entry_indentity<ET>(A, mid);
      combine_values(root, mid_et, false, op);
    }
    return Seq::node_join(P.first, P.second, root);
  }

  template <class BinaryOp, typename T>
  static node* multi_insert_sorted_count(ptr b, T* A, size_t n,
                                         BinaryOp const& op, size_t& sz) {
    // std::cout << "multi_insert_sorted" << std::endl;

    if (n == 0) return b.node_ptr();
    sz++;
    if (b.empty()) {
      std::cout << "empty " << n << std::endl;
      node* x = Seq::from_array(A, n);
      return x;
    }

    size_t tot = b.size() + n;
    // if (tot <= kBaseCaseSize) {
    if (b.is_compressed()) {
      // sz += count_size_ptr(std::move(b));
      // sz--;
      // sz++;
      return nullptr;
      // return multiinsert_bc(std::move(b), A, n, op);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    // node* rp = ptr::strip_flag(b.p);
    // auto root = Seq::cast_to_regular(rp);
    // ET e = root->entry.first;
    // ptr lc = ptr(root->lc, false);
    // ptr rc = ptr(root->rc, false);

    if (!root) root = Seq::single(e);

    K bk = Entry::get_key(e);
    // auto less_val = [&](ET& a) -> bool {
    //   return Entry::comp(Entry::get_key(a), bk);
    // };
    auto less_val = [&](auto const& a) -> bool {
      auto et = basic_node_helpers::get_entry_indentity<ET>(a);
      return Entry::comp(Entry::get_key(et), bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    // size_t mid = n / 2;
    // std::cout << e << " " << mid << " " << n << std::endl;
    // bool dup = (mid < n) && (!Entry::comp(bk, Entry::get_key(A[mid])));
    bool dup =
        (mid < n) &&
        (!Entry::comp(
            bk, Entry::get_key(
                    basic_node_helpers::get_entry_indentity<ET>(A, mid))));

    auto P = utils::fork<node*>(
        // true,  // Seq::do_parallel(b.size(), n),
        // !(mid == 0 || mid == n),
        0,
        // n >= 512,
        // Seq::do_parallel(b.size(), n),
        [&]() {
          return multi_insert_sorted_count(std::move(lc), A, mid, op, sz);
        },
        [&]() {
          return multi_insert_sorted_count(std::move(rc), A + mid + dup,
                                           n - mid - dup, op, sz);
        });

    // if (dup) combine_values(root, A[mid], false, op);
    if (dup) {
      auto mid_et = basic_node_helpers::get_entry_indentity<ET>(A, mid);
      combine_values(root, mid_et, false, op);
    }
    return Seq::node_join(P.first, P.second, root);
  }

  template <class BinaryOp, typename T>
  static node* multiinsert_bc_simple(node* b1, T* A, size_t n,
                                     BinaryOp const& op) {
    ET stack[Seq::size(b1) + 1];
    size_t offset = 0;
    auto copy_f = [&](ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(b1, copy_f);
    // Seq::decrement_recursive(n_b1);

    ET output[n + Seq::size(b1) + 1];

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      // auto const& k_b = Entry::get_key(A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      auto const& k_b = Entry::get_key(et);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        // parlay::move_uninitialized(output[out_off++], A[j]);
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        parlay::move_uninitialized(output[out_off++], et);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        // Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(A[j])));
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(et)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      // parlay::move_uninitialized(output[out_off++], A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      parlay::move_uninitialized(output[out_off++], et);
      j++;
    }

    assert(out_off == n + Seq::size(b1));

    // build returned tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      return Seq::make_compressed(output, out_off);
    }
  }

  template <class BinaryOp, typename T>
  static node* multi_insert_sorted_simple(node* b, T* A, size_t n,
                                          BinaryOp const& op, size_t& sz) {
    if (n == 0) return b;
    // sz++;

    if (Seq::size(b) == 0) {
      node* x = Seq::from_array(A, n);
      return x;
    }

    size_t tot = Seq::size(b) + n;
    if (tot <= kBaseCaseSize) {
      // return nullptr;
      return multiinsert_bc_simple(b, A, n, op);
    }

    // if (Seq::is_compressed(b)) {
    //   std::cout << Seq::size(b) << std::endl;
    // }
    // assert(!Seq::is_compressed(b));
    // auto [lc, e, rc, root] = Seq::expose(std::move(b));
    auto rt = Seq::is_regular(b) ? Seq::cast_to_regular(b)
                                 : Seq::to_tree((compressed_node*)b);
    // if (!rt) rt = Seq::single(e);
    assert(rt != nullptr);

    K bk = Entry::get_key(Seq::get_entry(rt));
    // auto less_val = [&](ET& a) -> bool {
    //   return Entry::comp(Entry::get_key(a), bk);
    // };
    auto less_val = [&](auto const& a) -> bool {
      auto et = basic_node_helpers::get_entry_indentity<ET>(a);
      return Entry::comp(Entry::get_key(et), bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    // std::cout << e << " " << mid << " " << n << std::endl;
    // bool dup = (mid < n) && (!Entry::comp(bk, Entry::get_key(A[mid])));
    bool dup =
        (mid < n) &&
        (!Entry::comp(
            bk, Entry::get_key(
                    basic_node_helpers::get_entry_indentity<ET>(A, mid))));

    auto P = utils::fork<node*>(
        // true,  // Seq::do_parallel(b.size(), n),
        !(mid == 0 || mid == n),
        // false,
        // Seq::do_parallel(b.size(), n),
        [&]() { return multi_insert_sorted_simple(rt->lc, A, mid, op, sz); },
        [&]() {
          return multi_insert_sorted_simple(rt->rc, A + mid + dup,
                                            n - mid - dup, op, sz);
        });

    // if (dup) combine_values(root, A[mid], false, op);
    if (dup) {
      auto mid_et = basic_node_helpers::get_entry_indentity<ET>(A, mid);
      combine_values(rt, mid_et, false, op);
    }
    return Seq::node_join(P.first, P.second, rt);
  }

  template <class VE, class CombineOp, class MapOp>
  static node* multiinsert_map_bc(ptr b1, std::pair<K, VE>* A, size_t n,
                                  CombineOp const& combine_op,
                                  MapOp const& map_op) {
    auto n_b1 = b1.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) { parlay::move_uninitialized(stack[offset++], a); };
    Seq::iterate_seq(n_b1, copy_f);
    assert(offset <= kBaseCaseSize);

    Seq::decrement_recursive(n_b1);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = A[j].first;
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        auto e = ET(A[j].first, map_op(A[j].second));
        parlay::move_uninitialized(output[out_off++], e);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, combine_op(Entry::get_val(re), A[j].second));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      auto e = ET(A[j].first, map_op(A[j].second));
      parlay::move_uninitialized(output[out_off++], e);
      j++;
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      auto ret = Seq::make_compressed(output, out_off);
      return ret;
    }
  }

  // Assumes the input is sorted and there are no duplicate keys
  template <class VE, class MapOp>
  static node* from_array_map(std::pair<K, VE>* A, size_t n,
                              MapOp const& map_op) {
    if (n <= 0) return Seq::empty();
    if (n == 1) return Seq::single(ET(A[0].first, map_op(A[0].second)));
    if (n >= B && n <= 2 * B) {
      ET stk[2 * B];
      parlay::parallel_for(
          0, n, [&](size_t i) { stk[i] = ET(A[i].first, map_op(A[i].second)); },
          1);
      return Seq::make_compressed(stk, n);
    }
    size_t mid = n / 2;

    regular_node* m =
        Seq::make_regular_node(ET(A[mid].first, map_op(A[mid].second)));

    auto P = utils::fork<node*>(
        n >= Seq::kNodeLimit, [&]() { return from_array_map(A, mid, map_op); },
        [&]() { return from_array_map(A + mid + 1, n - mid - 1, map_op); });

    return Seq::node_join(P.first, P.second, m);
  }

  // assumes array A is of length n and is sorted with no duplicates
  // combine: <V, VE> -> V
  // map: VE -> V
  template <class VE, class CombineOp, class MapOp>
  static node* multi_insert_sorted_map(ptr b, std::pair<K, VE>* A, size_t n,
                                       CombineOp const& combine_op,
                                       MapOp const& map_op) {
    if (b.empty()) {
      node* x = from_array_map(A, n, map_op);
      return x;
    }
    if (n == 0) return b.node_ptr();

    size_t tot = b.size() + n;
    // if (b.is_compressed() && tot <= kBaseCaseSize) {
    //   return multiinsert_map_bc(std::move(b), A, n, combine_op, map_op);
    // }
    if (tot <= Seq::kBaseCaseSize) {
      // if (tot <= Seq::kCompressionBlockSize) {
      return multiinsert_map_bc(std::move(b), A, n, combine_op, map_op);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    if (!root) root = Seq::single(e);

    K bk = Entry::get_key(e);
    auto less_val = [&](std::pair<K, VE>& a) -> bool {
      return Entry::comp(a.first, bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid].first));

    auto P = utils::fork<node*>(
        true,  // Seq::do_parallel(b.size(), n),
        [&]() {
          return multi_insert_sorted_map(std::move(lc), A, mid, combine_op,
                                         map_op);
        },
        [&]() {
          return multi_insert_sorted_map(std::move(rc), A + mid + dup,
                                         n - mid - dup, combine_op, map_op);
        });

    if (dup) combine_values_v(root, A[mid].second, combine_op);
    return Seq::node_join(P.first, P.second, root);
  }

  template <class VE, class BinaryOp>
  static bool multiupdate_sorted_inplace_bc(node* b, std::pair<K, VE>* A,
                                            size_t n, BinaryOp const& op) {
    size_t ct = 0;
    auto update_f = [&](ET et) {
      K key_et = Entry::get_key(et);
      V ret_v = Entry::get_val(et);
      if (ct < n) {
        K key_a = A[ct].first;
        while (Entry::comp(key_a, key_et) && ct < n - 1) {
          ct++;
          key_a = A[ct].first;
        }
        if (!Entry::comp(key_a, key_et) && !Entry::comp(key_et, key_a)) {
          ret_v = op(Entry::get_val(et), A[ct].second);
          ct++;
        }
      }
      return ret_v;
    };
    Seq::inplace_update(b, update_f);
    return true;
  }

  template <class VE, class BinaryOp>
  static bool multi_update_sorted_inplace(node* b, std::pair<K, VE>* A,
                                          size_t n, BinaryOp const& op) {
    if (b == nullptr) return true;
    if (n == 0) return true;
    assert(Seq::ref_cnt(b) == 1);

    if (Seq::is_compressed(b)) {
      return multiupdate_sorted_inplace_bc(b, A, n, op);
    }

    auto rb = (regular_node*)b;

    K bk = get_key(b);
    auto less_val = [&](std::pair<K, VE>& a) -> bool {
      return Entry::comp(a.first, bk);
    };

    size_t mid = utils::PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid].first));

    auto P = utils::fork<bool>(
        true,
        [&]() {
          return multi_update_sorted_inplace(std::move(rb->lc), A, mid, op);
        },
        [&]() {
          return multi_update_sorted_inplace(std::move(rb->rc), A + mid + dup,
                                             n - mid - dup, op);
        });

    if (dup) combine_values_v(rb, A[mid].second, op);
    return P.first && P.second;
  }

  static bool multifind_sorted_bc(node* b, K* A, size_t n, V* ret,
                                  size_t offset) {
    size_t ct = 0;
    auto emit_f = [&](const ET& et) {
      K key_et = Entry::get_key(et);
      if (ct < n) {
        K key_a = A[ct];
        while (Entry::comp(key_a, key_et) && (ct < n - 1)) {
          ct++;
          key_a = A[ct];
        }
        if (!Entry::comp(key_a, key_et) && !Entry::comp(key_et, key_a)) {
          ret[offset + ct] = Entry::get_val(et);
          ct++;
        }
      }
    };
    Seq::iterate_seq(b, emit_f);
    return true;
  }

  static bool multi_find_sorted(node* b, K* A, size_t n, V* ret,
                                size_t offset) {
    if (!b) return true;
    if (n == 0) return true;

    if (Seq::is_compressed(b)) {
      return multifind_sorted_bc(b, A, n, ret, offset);
    }

    auto rb = (regular_node*)b;
    K bk = get_key(rb);
    auto less_val = [&](K a) -> bool { return Entry::comp(a, bk); };
    size_t mid =
        parlay::internal::binary_search(parlay::make_slice(A, A + n), less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));
    utils::fork<bool>(
        true, [&]() { return multi_find_sorted(rb->lc, A, mid, ret, offset); },
        [&]() {
          return multi_find_sorted(rb->rc, A + mid + dup, n - mid - dup, ret,
                                   offset + mid + dup);
        });

    if (dup) ret[offset + mid] = get_val(rb);
    return true;
  }

  template <class InTree, class Func>
  static node* map(typename InTree::ptr b, Func const& f) {
    auto g = [&](typename InTree::ET& a) {
      return ET(InTree::Entry::get_key(a), f(a));
    };
    return Seq::template map<InTree>(std::move(b), g);
  }

  template <class InTree, class Func>
  static node* map_set(typename InTree::ptr b, Func const& f) {
    return Seq::template map<InTree>(std::move(b), f);
  }

  // TODO
  //  template<class Seq1, class Func>
  //  static node* map_filter(typename Seq1::node* b, const Func& f,
  //                          size_t granularity=utils::node_limit) {
  //    auto g = [&] (typename Seq1::ET& a) {
  //      std::optional<V> v = f(a);
  //      if (v) return std::optional<ET>(ET(Seq1::Entry::get_key(a),*v));
  //      else return std::optional<ET>();
  //    };
  //
  //    return Seq::template map_filter<Seq1>(b, g, granularity);
  //  }

  template <class BinaryOp, typename T>
  static node* multiinsert_bc_append(ptr b1, T* A, size_t n,
                                     BinaryOp const& op) {
    assert(b1.size() + n <= 2 * B);
    auto n_b1 = b1.node_ptr();
    size_t offset = Seq::size(n_b1);
    auto c = Seq::cast_to_compressed(n_b1);
    uint8_t* data_start =
        (((uint8_t*)c) + sizeof(typename Seq::aug_compressed_node));
    ET* stack = (ET*)data_start;
    for (size_t i = 0; i < n; i++) {
      ET et = basic_node_helpers::get_entry_indentity<ET>(A, i);
      stack[offset + i] = et;
      c->aug_val = Entry::combine(c->aug_val, Entry::from_entry(et));
    }
    c->s = offset + n;
    return c;
  }

  template <class BinaryOp, typename T>
  static node* multiinsert_bc_inplace(ptr b1, T* A, size_t n,
                                      BinaryOp const& op) {
    auto n_b1 = b1.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    Seq::decrement_recursive(n_b1);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      // auto const& k_b = Entry::get_key(A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      auto const& k_b = Entry::get_key(et);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        // parlay::move_uninitialized(output[out_off++], A[j]);
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        parlay::move_uninitialized(output[out_off++], et);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        // Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(A[j])));
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(et)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      // parlay::move_uninitialized(output[out_off++], A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      parlay::move_uninitialized(output[out_off++], et);
      j++;
    }

    // build returned tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      return Seq::make_compressed(output, out_off);
    }
  }

  // template <class BinaryOp>
  // static node* multiinsert_bc(ptr b1, ET* A, size_t n, BinaryOp const& op)
  // {
  template <class BinaryOp, typename T>
  static node* multiinsert_bc(ptr b1, T* A, size_t n, BinaryOp const& op) {
    auto n_b1 = b1.node_ptr();

    ET stack[Seq::size(n_b1) + 1];
    // auto stack = parlay::sequence<ET>::uninitialized(Seq::size(n_b1) + 1);
    size_t offset = 0;
    auto copy_f = [&](ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    Seq::decrement_recursive(n_b1);

    // auto output = parlay::sequence<ET>::uninitialized(offset + n + 1);
    ET output[offset + n + 1];
    // size_t out_off = offset;

    // merge
    size_t nA = offset;
    size_t nB = n;
    size_t i = 0, j = 0, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      // auto const& k_b = Entry::get_key(A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      auto const& k_b = Entry::get_key(et);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        // parlay::move_uninitialized(output[out_off++], A[j]);
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        parlay::move_uninitialized(output[out_off++], et);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        // Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(A[j])));
        auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
        Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(et)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      // parlay::move_uninitialized(output[out_off++], A[j]);
      auto et = basic_node_helpers::get_entry_indentity<ET>(A, j);
      parlay::move_uninitialized(output[out_off++], et);
      j++;
    }

    // build returned tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
      // return Seq::to_tree_impl(output.data(), out_off);
    } else {
      // return Seq::make_compressed(output, out_off);
      return Seq::from_array(output, out_off);
    }
  }

  template <class BinaryOp>
  static node* union_bc(ptr b1, ptr b2, BinaryOp const& op) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) {  // TODO: copy or ref?
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq::iterate_seq(n_b2, copy_f);
    assert(offset <= kBaseCaseSize);

    Seq::decrement_recursive(n_b1);
    Seq::decrement_recursive(n_b2);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2;
    size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        parlay::move_uninitialized(output[out_off++], stack[j]);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, op(Entry::get_val(stack[j]), Entry::get_val(re)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      parlay::move_uninitialized(output[out_off++], stack[j]);
      j++;
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  static node* difference_bc(ptr b1, ptr b2) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) {  // TODO: copy or ref?
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq::iterate_seq(n_b2, copy_f);
    assert(offset <= kBaseCaseSize);

    Seq::decrement_recursive(n_b1);
    Seq::decrement_recursive(n_b2);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2;
    size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        j++;
      } else {
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  template <class Seq1, class Seq2, class BinaryOp>
  static node* intersect_bc(typename Seq1::ptr b1, typename Seq2::ptr b2,
                            BinaryOp const& op) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&](ET a) {  // TODO: copy or ref?
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq1::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq2::iterate_seq(n_b2, copy_f);
    assert(offset <= kBaseCaseSize);

    Seq1::decrement_recursive(n_b1);
    Seq2::decrement_recursive(n_b2);

    ET output[kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2;
    size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      auto const& k_a = Entry::get_key(stack[i]);
      auto const& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        i++;
      } else if (comp(k_b, k_a)) {
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, op(Entry::get_val(stack[j]), Entry::get_val(re)));
        out_off++;
        i++;
        j++;
      }
    }

    // build tree
    if (out_off < B) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }
};

}  // namespace cpam
