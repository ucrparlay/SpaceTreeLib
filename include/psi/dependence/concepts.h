#ifndef PSI_DEPENDENCE_CONCEPTS_H_
#define PSI_DEPENDENCE_CONCEPTS_H_

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <variant>

namespace psi
{

// NOTE: in memory of SFINE
// template<typename>
// struct is_pair : std::false_type {};
// template<typename T, typename U>
// struct is_pair<std::pair<T, U>> : std::true_type {};
// template<typename T>
// concept is_pair = is_pair<T>::value;

template <typename T, typename node>
concept check_pointer_to_node =
	std::is_pointer_v<T> &&
	std::is_base_of_v<node, std::remove_pointer_t<T>>;

template <typename T>
concept is_pair = requires {
	requires std::same_as<
		T, std::pair<typename T::first_type, typename T::second_type>>;
};

template <typename T, typename Point>
concept is_box = requires {
	requires is_pair<T> &&
			 std::same_as<typename T::first_type,
				      typename Point::bp_type> &&
			 std::same_as<typename T::second_type,
				      typename Point::bp_type>;
};

// NOTE: tag whether a node has non-trivial augmentation
template <typename T>
concept node_has_non_trivial_aug =
	(!std::is_empty_v<T>) && (!std::same_as<T, std::monostate>);

// NOTE:: Helper function to find and return the variable of the specified type
template <typename T, typename First, typename... Rest>
constexpr T &find_var(First &first, Rest &...rest)
{
	if constexpr (std::same_as<T, First>) {
		return first;
	} else {
		return find_var<T>(rest...);
	}
}

// NOTE: Overload for the case when the type is not found
template <typename T>
constexpr T &find_var()
{
	static_assert(!std::same_as<T, T>, "Type not found in parameter pack");
	return *reinterpret_cast<T *>(nullptr);
}

// Example usage
// template<typename T, typename... Args>
//     requires ContainsType<T, Args...>
// consteval T& getVariable(Args&... args) {
//     return findVariable<T>(args...);
// }

// NOTE: Concept to check if a function can be called with a parameter
template <typename ReturnType, typename Func, typename Arg>
concept callable_with_arg = requires(Func func, Arg arg) {
	{ func(arg) } -> std::convertible_to<ReturnType>;
};

// NOTE: Function to handle the conditional call
// if func accepts a parameter arg, then call it; otherwise return func()
// directly
template <typename ReturnType, typename Func, typename Arg>
ReturnType invoke_with_optional_arg(Func &&func, Arg &&arg)
{
	if constexpr (callable_with_arg<ReturnType, Func, Arg>) {
		return func(std::forward<Arg>(arg));
	} else {
		return func();
	}
}

// NOTE: define the what is a binary node
template <typename T,
	  template <typename, typename, typename> typename BinaryNodeT>
concept check_binary_node =
	std::is_base_of_v<BinaryNodeT<typename T::pt_type, typename T::st_type,
				      typename T::at_type>,
			  T>;

// NOTE: helper for decide a multi-way node
template <typename T,
	  template <typename, uint_fast8_t, typename,
		    typename> typename MultiNodeT,
	  uint_fast8_t... Ns>
concept is_multi_node_helper =
	(std::is_base_of_v<MultiNodeT<typename T::pt_type, Ns,
				      typename T::st_type, typename T::at_type>,
			   T> ||
	 ...);

// NOTE: define the what is a multi-way node
template <typename T, template <typename, uint_fast8_t, typename,
				typename> typename MultiNodeT>
concept check_multi_node =
	is_multi_node_helper<T, MultiNodeT, 2, 3, 4, 5, 6, 7, 8, 9, 10>;

// NOTE: define the what is a dynamic node
template <typename T,
	  template <typename, typename, typename> typename DynamicNodeT>
concept check_dynamic_node =
	std::is_base_of_v<DynamicNodeT<typename T::pt_type, typename T::st_type,
				       typename T::at_type>,
			  T>;

// NOTE: tag aug point
template <typename T>
concept is_aug_point = requires(T t) {
	{ t.aug_point_tag() } -> std::same_as<void>;
};

/*
 * Carries a bounding box. Tests get_box(), which is what the library calls; a
 * member named box is not enough, and used to be the whole test, so an
 * augmentation that renamed it silently fell back to the hyperplane paths.
 */
template <typename T>
concept has_box = requires(T const &t) {
	{ t.get_box() };
};

// NOTE: tag a orth tree
template <typename T>
concept is_orth_tree = requires(T t) {
	{ t.orth_tree_tag() } -> std::same_as<void>;
};

// NOTE: tag a kd tree
template <typename T>
concept is_kd_tree = requires(T t) {
	{ t.kd_tree_tag() } -> std::same_as<void>;
};

template <typename T>
concept is_p_tree = requires(T t) {
	{ t.p_tree_tag() } -> std::same_as<void>;
};

template <typename T>
concept is_object_median_split = requires(T t) {
	{ t.object_median_tag() } -> std::same_as<void>;
};

template <typename T>
concept is_spatial_median_split = requires(T t) {
	{ t.spatial_median_tag() } -> std::same_as<void>;
};

/*
 * What kd_tree and orth_tree require of a user supplied augmentation. create
 * and update are not here: they take the very interior_type that holds the
 * augmentation, which is not nameable at the point these are checked.
 */
template <typename A, typename slice_type>
concept leaf_augmentation = requires(A a, slice_type in) {
	A();
	A(in);
	a.update_aug(in);
	a.reset();
};

template <typename A>
concept interior_augmentation = requires(A a, bool flag, size_t n) {
	a.set_parallel_flag(flag);
	a.reset_parallel_flag();
	{ a.get_parallel_flag_ini_status() } -> std::convertible_to<bool>;
	{ a.force_parallel(n) } -> std::convertible_to<bool>;
};

template <typename T>
concept supports_force_parallel = requires(T t) {
	{ t.force_parallel() } -> std::same_as<bool>;
};

template <typename T>
concept is_max_stretch_dim = requires(T t) {
	{ t.max_stretch_tag() } -> std::same_as<void>;
};

template <typename T>
concept is_rotate_dim_split = requires(T t) {
	{ t.rotate_dim_tag() } -> std::same_as<void>;
};

} // namespace psi

#endif // PSI_DEPENDENCE_CONCEPTS_H_
