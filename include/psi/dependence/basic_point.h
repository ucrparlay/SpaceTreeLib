#ifndef PSI_DEPENDENCE_BASIC_POINT_H_
#define PSI_DEPENDENCE_BASIC_POINT_H_

#include <cstdint>
#include <tuple>
#include <type_traits>
#include <variant>

#include "comparator.h"
#include "parlay/alloc.h"
#include "parlay/delayed.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/type_traits.h"
#include "parlay/utilities.h"

namespace psi
{

template <typename T, uint_fast8_t d>
struct basic_point {
	using coord_type = T;
	using coords_type = std::array<T, d>;
	using dis_type = std::conditional_t<std::is_integral_v<coord_type>,
					    int_fast64_t, double>;
	using num_type = num_comparator<coord_type>;
	using dims_type = uint_fast8_t;
	using bp_type = basic_point<T, d>;

	basic_point()
	{
	}

	explicit basic_point(T const val)
	{
		this->pnt.fill(val);
	}

	explicit basic_point(coords_type const &_pnt) : pnt(_pnt)
	{
	}

	explicit basic_point(parlay::slice<T *, T *> x)
	{
		assert(x.size() == d);
		for (dims_type i = 0; i < d; ++i) {
			this->pnt[i] = x[i];
		}
	}

	explicit basic_point(T const *x)
	{
		for (dims_type i = 0; i < d; ++i) {
			this->pnt[i] = x[i];
		}
	}

	basic_point const min_coords(basic_point const &b) const
	{
		basic_point p;

		if constexpr (d == 2) {
			p.pnt[0] = num_type::min(this->pnt[0], b.pnt[0]);
			p.pnt[1] = num_type::min(this->pnt[1], b.pnt[1]);
		} else if constexpr (d == 3) {
			p.pnt[0] = num_type::min(this->pnt[0], b.pnt[0]);
			p.pnt[1] = num_type::min(this->pnt[1], b.pnt[1]);
			p.pnt[2] = num_type::min(this->pnt[2], b.pnt[2]);
		} else {
			for (dims_type i = 0; i < d; ++i) {
				p.pnt[i] =
					num_type::min(this->pnt[i], b.pnt[i]);
			}
		}

		return p;
	}

	basic_point const max_coords(basic_point const &b) const
	{
		basic_point p;

		if constexpr (d == 2) {
			p.pnt[0] = num_type::max(this->pnt[0], b.pnt[0]);
			p.pnt[1] = num_type::max(this->pnt[1], b.pnt[1]);
		} else if constexpr (d == 3) {
			p.pnt[0] = num_type::max(this->pnt[0], b.pnt[0]);
			p.pnt[1] = num_type::max(this->pnt[1], b.pnt[1]);
			p.pnt[2] = num_type::max(this->pnt[2], b.pnt[2]);
		} else {
			for (dims_type i = 0; i < d; ++i) {
				p.pnt[i] =
					num_type::max(this->pnt[i], b.pnt[i]);
			}
		}

		return p;
	}

	static consteval auto get_dim()
	{
		return d;
	}

	coords_type &get_coords()
	{
		return pnt;
	}
	coords_type const &get_coords() const
	{
		return pnt;
	}

	bool same_dimension(basic_point const &b) const
	{
		return *this == b;
	}

	bool operator==(basic_point const &x) const
	{
		for (dims_type i = 0; i < d; ++i) {
			if (!num_type::eq(this->pnt[i], x.pnt[i]))
				return false;
		}
		return true;
	}

	bool operator<(basic_point const &x) const
	{
		for (dims_type i = 0; i < d; ++i) {
			if (num_type::lt(this->pnt[i], x.pnt[i]))
				return true;
			else if (num_type::gt(this->pnt[i], x.pnt[i]))
				return false;
			else
				continue;
		}
		return false;
	}

	coord_type &operator[](dims_type i)
	{
		return pnt[i];
	}

	coord_type const &operator[](dims_type i) const
	{
		return pnt[i];
	}

	friend std::ostream &operator<<(std::ostream &o, basic_point const &a)
	{
		o << "(";
		for (dims_type i = 0; i < d; ++i) {
			o << a.pnt[i] << (i == d - 1 ? "" : ", ");
		}
		o << ") " << std::flush;
		return o;
	}

	coords_type pnt;
};

template <typename T, uint_fast8_t d, class AugType = std::monostate>
	requires parlay::is_less_than_comparable_v<AugType>
struct aug_point : basic_point<T, d> {
	using bp_type = basic_point<T, d>;
	using dims_type = bp_type::dims_type;
	using at_type = AugType;

	static consteval bool is_non_trivial_augmentation()
	{
		return !std::is_same_v<AugType, std::monostate>;
	}

	void aug_point_tag()
	{
	}

	aug_point()
	{
	}

	explicit aug_point(T const val) : basic_point<T, d>(val)
	{
	}

	explicit aug_point(typename basic_point<T, d>::coords_type const &_pnt)
	    : basic_point<T, d>(_pnt)
	{
	}

	explicit aug_point(parlay::slice<T *, T *> x) : basic_point<T, d>(x)
	{
	}

	explicit aug_point(T const *x) : basic_point<T, d>(x)
	{
	}

	aug_point(typename basic_point<T, d>::coords_type const &_pnt,
		  AugType const &aug)
	    : basic_point<T, d>(_pnt), aug(aug)
	{
	}

	aug_point(parlay::slice<T *, T *> x, AugType const &aug)
	    : basic_point<T, d>(x), aug(aug)
	{
	}

	aug_point(T const *x, AugType const &aug)
	    : basic_point<T, d>(x), aug(aug)
	{
	}

	aug_point(aug_point const &p) : basic_point<T, d>(p), aug(p.aug)
	{
	}

	aug_point(aug_point &&p) noexcept
	    : basic_point<T, d>(std::move(p)), aug(std::move(p.aug))
	{
	}

	aug_point &operator=(aug_point const &p)
	{
		basic_point<T, d>::operator=(p);
		aug = p.aug;
		return *this;
	}

	aug_point &operator=(aug_point &&p) noexcept
	{
		basic_point<T, d>::operator=(std::move(p));
		aug = std::move(p.aug);
		return *this;
	}

	bool operator==(aug_point const &x) const
	{
		return basic_point<T, d>::operator==(x) && aug == x.aug;
	}

	bool operator<(aug_point const &x) const
	{
		return basic_point<T, d>::operator<(x)
			       ? true
			       : (basic_point<T, d>::operator==(x) ? aug < x.aug
								   : false);
	}

	friend std::ostream &operator<<(std::ostream &o, aug_point const &a)
	{
		o << "[" << static_cast<basic_point<T, d> const &>(a) << " -> "
		  << a.aug << "] " << std::flush;
		return o;
	}

	bp_type min_coords(aug_point const &b) const
	{
		return static_cast<bp_type const &>(*this).min_coords(
			static_cast<bp_type const &>(b));
	}

	bp_type max_coords(aug_point const &b) const
	{
		return static_cast<bp_type const &>(*this).max_coords(
			static_cast<bp_type const &>(b));
	}

	AugType &get_aug()
	{
		return aug;
	}

	AugType const &get_aug() const
	{
		return aug;
	}

	template <typename... Args>
	void set_aug_member(Args &&...args)
	{
		return aug.set_member(std::forward<Args>(args)...);
	}

	AugType aug;
};

} // namespace psi

#endif // PSI_DEPENDENCE_BASIC_POINT_H_
