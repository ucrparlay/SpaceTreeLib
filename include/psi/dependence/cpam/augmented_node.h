#pragma once
#include "basic_node.h"

// *******************************************
//   AUGMENTED NODE
// *******************************************
namespace psi
{
namespace cpam
{

// Creates an augmented node from a basic node.
// The new augmented entry is a pair of the original entry and the aumented
// value.   The Entry struct must have the static inteface:
//   entry_t;
//   aug_t;
//   get_empty() -> aug_t;
//   from_entry(entry_t) -> aug_t;
//   combine(aut_t, aug_t) -> aug_t;
template <class balance, class Entry, class AugEntryEncoder, size_t block_size>
struct aug_node
    : public basic_node<
	      balance,
	      std::pair<typename Entry::entry_t, typename Entry::aug_t>,
	      null_encoder, block_size> {
	using at_type = typename Entry::aug_t;
	using et_type = typename Entry::entry_t;
	using aug = aug_node<balance, Entry, AugEntryEncoder, block_size>;
	using basic = basic_node<balance, std::pair<et_type, at_type>,
				 null_encoder, block_size>;
	using node = typename basic::node;
	using regular_node = typename basic::regular_node;
	using compressed_node = typename basic::compressed_node;
	using allocator = typename basic::allocator;
	using complex_allocator = typename basic::complex_allocator;

	using basic::B;
	using basic::check_compressed_node;
	using basic::empty;
	using basic::increment_count;
	using basic::node_limit;
	using basic::ref_cnt;
	using basic::size;
	using basic::will_be_compressed;

	static et_type &get_entry(node *a)
	{
		return basic::cast_to_regular(a)->entry.first;
	}
	static et_type *get_entry_p(node *a)
	{
		return &basic::cast_to_regular(a)->entry.first;
	}
	static void set_entry(node *a, et_type const &e)
	{
		basic::cast_to_regular(a)->entry.first = e;
	}

	struct aug_compressed_node {
		node_size_t r; // reference count, top-bit is always "0"
		node_size_t s; // number of entries used (size)
		node_size_t size_in_bytes;
		at_type aug_val;
		bool is_sorted = true;
	};

	static bool is_regular(node *a)
	{
		return !a || ((regular_node *)a)->r & basic::top_bit;
	}
	static bool is_compressed(node *a)
	{
		return !is_regular(a);
	}
	static regular_node *cast_to_regular(node *a)
	{
		assert(is_regular(a));
		return (regular_node *)a;
	}
	static aug_compressed_node *cast_to_compressed(node *a)
	{
		assert(is_compressed(a));
		return (aug_compressed_node *)a;
	}

	// TODO: include case for compressed
	static at_type aug_val(node *a)
	{
		if (a == NULL)
			return Entry::get_empty();
		if (basic::is_regular(a))
			return (basic::cast_to_regular(a)->entry).second;
		auto c = (aug_compressed_node *)a;
		return c->aug_val;
	}

	static at_type const &aug_val_ref(node *a)
	{
		if (a == NULL) {
			/* get_empty() returns by value; a reference bound to
			 * it dies with the return statement. */
			static at_type const empty = Entry::get_empty();
			return empty;
		}
		if (basic::is_regular(a))
			return (basic::cast_to_regular(a)->entry).second;
		auto c = (aug_compressed_node *)a;
		return c->aug_val;
	}

	// updates augmented value using f, instead of recomputing
	template <class F>
	static void lazy_update(node *ov, F f)
	{
		regular_node *a = basic::cast_to_regular(ov);
		basic::update(a);
		(a->entry).second = f((a->entry).second);
	}

	static void update(node *ov)
	{
		regular_node *a = basic::cast_to_regular(ov);
		basic::update(a);
		at_type av = Entry::from_entry(get_entry(a));
		if (a->lc)
			av = Entry::combine(aug_val(a->lc), std::move(av));
		if (a->rc)
			av = Entry::combine(std::move(av), aug_val(a->rc));
		(a->entry).second = std::move(av);
	}

	// Handles both regular and compressed nodes.
	static void free_node(node *va)
	{
		if (basic::is_regular(va)) {
			auto a = basic::cast_to_regular(va);
			std::get<0>((a->entry)).~et_type();
			std::get<1>((a->entry)).~at_type();
			allocator::free(a);
		} else {
			auto c = cast_to_compressed(va);
			uint8_t *data_start =
				(((uint8_t *)c) + sizeof(aug_compressed_node));
			c->aug_val.~at_type();
			AugEntryEncoder::destroy(data_start, c->s);
			auto array_size = c->size_in_bytes;
			utils::free_array<uint8_t>((uint8_t *)va, array_size);
		}
	}

	// precondition: a is non-null, and caller has a reference count on a.
	static bool decrement_count(node *a)
	{
		if ((utils::fetch_and_add(&basic::generic_node(a)->r, -1) &
		     basic::low_bit_mask) == 1) {
			free_node(a);
			return true;
		}
		return false;
	}

	// atomically decrement ref count and deletes node if zero
	static bool decrement(node *t)
	{
		if (t) {
			return decrement_count(t);
		}
		return false;
	}

	// atomically decrement ref count and if zero:
	//   delete node and recursively decrement the two children
	static void decrement_recursive(node *t)
	{
		if (!t)
			return;
		if (basic::is_regular(t)) {
			node *lsub = basic::cast_to_regular(t)->lc;
			node *rsub = basic::cast_to_regular(t)->rc;
			if (decrement(t)) {
				utils::fork_no_result(
					basic::size(lsub) >= node_limit,
					[&]() { decrement_recursive(lsub); },
					[&]() { decrement_recursive(rsub); });
			}
		} else {
			decrement(t);
		}
	}

	template <typename T>
	// static regular_node* make_regular_node(et_type e) {
	static regular_node *make_regular_node(T const &e)
	{
		std::pair<et_type, at_type> ea;
		// ea.first = e;
		ea.first = basic_node_helpers::get_entry_indentity<et_type>(e);
		return basic::make_regular_node(ea);
	}

	static regular_node *single(et_type const &e)
	{
		at_type av = Entry::from_entry(e);
		return basic::single(std::make_pair(e, av));
	}

	// static regular_node* single(const et_type& e) {
	//   at_type av = Entry::from_entry(e);
	//   return basic::single(std::make_pair(e, av));
	// }

	template <typename F>
	static void inplace_update(node *a, F const &f)
	{
		assert(!is_regular(a));
		auto c = cast_to_compressed(a);
		uint8_t *data_start =
			(((uint8_t *)c) + sizeof(aug_compressed_node));
		AugEntryEncoder::inplace_update(data_start, c->s, f);
	}

	static void reorder(node *p)
	{
		auto c = cast_to_compressed(p);
		uint8_t *data_start =
			(((uint8_t *)c) + sizeof(aug_compressed_node));
		if (!c->is_sorted) {
			AugEntryEncoder::reorder(data_start, c->s);
			c->is_sorted = true;
		}
	}

	template <typename F, bool ReorderCompress = true>
	static void iterate_seq(node *a, F const &f)
	{
		if (!a)
			return;
		if (is_regular(a)) {
			auto r = cast_to_regular(a);
			iterate_seq(r->lc, f);
			f(get_entry(r));
			iterate_seq(r->rc, f);
		} else {
			auto c = cast_to_compressed(a);
			uint8_t *data_start =
				(((uint8_t *)c) + sizeof(aug_compressed_node));
			if constexpr (ReorderCompress) {
				reorder(c);
			}
			AugEntryEncoder::decode(data_start, c->s, f);
		}
	}

	template <typename F, bool ReorderCompress = true>
	static bool iterate_cond(node *a, F const &f)
	{
		if (!a)
			return true;
		if (is_regular(a)) {
			auto r = cast_to_regular(a);
			bool ret = iterate_cond(r->lc, f);
			if (!ret)
				return ret;
			ret = f(get_entry(r));
			if (!ret)
				return ret;
			ret = iterate_cond(r->rc, f);
			return ret;
		} else {
			auto c = cast_to_compressed(a);
			uint8_t *data_start =
				((uint8_t *)c) + sizeof(aug_compressed_node);
			if constexpr (ReorderCompress) {
				reorder(c);
			}
			return AugEntryEncoder::decode_cond(data_start, c->s,
							    f);
		}
	}

	template <class F, class Comp, class K>
	static std::optional<et_type>
	find_compressed(node *b, F const &f, Comp const &comp, K const &k)
	{
		auto c = cast_to_compressed(b);
		uint8_t *data_start =
			((uint8_t *)c) + sizeof(aug_compressed_node);
		return AugEntryEncoder::find(data_start, c->s, f, comp, k);
	}

	static node *finalize(node *root)
	{
		auto sz = basic::size(root);
		// std::cout << "finalzie in aug node" << std::endl;
		// std::cout << "sz " << sz << ", B " << B << std::endl;
		assert(sz > 0);
		if (sz < B) {
			auto ret = make_compressed_node(root);
			// std::cout << "finalize: " << sz << " " << B <<
			// std::endl;
			decrement_recursive(root);
			return ret;
		}
		return root;
	}

	// Used by gc_type to copy a compressed node. TODO: update to work
	// correctly with diff-encoding.
	static node *make_compressed_node(node *b)
	{
		return basic_node_helpers::make_compressed_node<aug>(b, B);
	}

	// takes a pointer to an array of ETs, and a length of the number of ETs
	// to construct, and returns a compressed node.
	template <typename T>
	// static compressed_node* make_single_compressed_node(et_type* e,
	// size_t s)
	// {
	static compressed_node *make_single_compressed_node(T *e, size_t s)
	{
		assert(s <= 2 * B);

		/* Sized for 2 * B entries so a later append fits in place,
		 * though only s are written below. default_entry_encoder
		 * ignores the pointer; a compressing one would need e. */
		size_t encoded_size = AugEntryEncoder::encoded_size(
			static_cast<et_type *>(nullptr), 2 * B);
		size_t node_size = sizeof(aug_compressed_node) + encoded_size;
		aug_compressed_node *c_node = (aug_compressed_node *)
			utils::new_array_no_init<uint8_t>(node_size);

		c_node->r = 1;
		c_node->s = s;
		c_node->size_in_bytes = node_size;
		c_node->is_sorted = true;

		uint8_t *encoded_data =
			(((uint8_t *)c_node) + sizeof(aug_compressed_node));
		parlay::assign_uninitialized(
			c_node->aug_val,
			AugEntryEncoder::encode(e, s, encoded_data));

		return (compressed_node *)c_node;
	}

	static node *make_compressed(node *l, node *r, regular_node *e)
	{
		return basic_node_helpers::make_compressed<aug>(l, r, e, B);
	}

	template <typename T>
	// static node* make_compressed(et_type* stack, size_t tot) {
	static node *make_compressed(T *stack, size_t tot)
	{
		return basic_node_helpers::make_compressed<aug>(stack, tot, B);
	}

	static node *make_compressed(node *l, node *r)
	{
		return make_compressed(l, r, nullptr);
	}

	static et_type *compressed_node_elms(void *ac, et_type *tmp_arr)
	{
		auto c = (aug_compressed_node *)ac;
		uint8_t *data_start =
			(((uint8_t *)c) + sizeof(aug_compressed_node));
		size_t i = 0;
		auto f = [&](et_type const &et) {
			parlay::assign_uninitialized(tmp_arr[i++], et);
		};
		AugEntryEncoder::decode(data_start, c->s, f);
		return tmp_arr;
	}

	// F applied to entries.
	template <class F>
	static size_t size_in_bytes(node *a, F const &f)
	{
		size_t total = 0;
		if (!a)
			return total;
		if (basic::is_compressed(a)) {
			auto c = basic::cast_to_compressed(a);
			total += c->size_in_bytes;
			auto fn = [&](auto const &et) { total += f(et); };
			iterate_seq(c, fn);
		} else {
			auto r = basic::cast_to_regular(a);
			total += size_in_bytes(r->lc, f);
			total += size_in_bytes(r->rc, f);
			total += sizeof(regular_node);
			total += f(get_entry(r));
		}
		return total;
	}
};

} // namespace cpam
} // namespace psi
