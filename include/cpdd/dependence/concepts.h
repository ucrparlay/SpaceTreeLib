
#pragma once
#include <concepts>
#include <tuple>
#include <type_traits>

#include "tree_node.h"

namespace cpdd {

// NOTE: In memory of SFINE
// template<typename>
// struct is_pair : std::false_type {};
// template<typename T, typename U>
// struct is_pair<std::pair<T, U>> : std::true_type {};
// template<typename T>
// concept IsPair = is_pair<T>::value;

template <typename T>
concept IsPair = requires {
  requires std::same_as<
      T, std::pair<typename T::first_type, typename T::second_type>>;
};

template <typename T, typename Point>
concept IsBox = requires {
  requires IsPair<T> && std::same_as<typename T::first_type, Point> &&
               std::same_as<typename T::second_type, Point>;
};

template <typename T>
concept IsPointer = std::is_pointer_v<T>;

// NOTE:  Concept to check if a type is present in a parameter pack
template <typename T, typename... Args>
concept ContainsType = (std::same_as<T, Args> || ...);

// NOTE:: Helper function to find and return the variable of the specified type
template <typename T, typename First, typename... Rest>
constexpr T& FindVar(First& first, Rest&... rest) {
  if constexpr (std::same_as<T, First>) {
    return first;
  } else {
    return FindVar<T>(rest...);
  }
}

// NOTE: Overload for the case when the type is not found
template <typename T>
constexpr T& FindVar() {
  static_assert(!std::same_as<T, T>, "Type not found in parameter pack");
  return *reinterpret_cast<T*>(nullptr);
}

// Example usage
// template<typename T, typename... Args>
//     requires ContainsType<T, Args...>
// consteval T& getVariable(Args&... args) {
//     return findVariable<T>(args...);
// }

// NOTE: Concept to check if a function can be called with a parameter
template <typename ReturnType, typename Func, typename Arg>
concept CallableWithArg = requires(Func func, Arg arg) {
  { func(arg) } -> std::convertible_to<ReturnType>;
};

// NOTE: Function to handle the conditional call
template <typename ReturnType, typename Func, typename Arg>
ReturnType InvokeWithOptionalArg(Func&& func, Arg&& arg) {
  if constexpr (CallableWithArg<ReturnType, Func, Arg>) {
    return func(std::forward<Arg>(arg));
  } else {
    return func();
  }
}

template <typename T>
concept IsBinaryNode = std::is_base_of_v<
    BinaryNode<typename T::PT, typename T::ST, typename T::AT>, T>;

template <typename T>
concept IsMultiNode =
    std::is_base_of_v<
        MultiNode<typename T::PT, 2, typename T::ST, typename T::AT>, T> ||
    std::is_base_of_v<
        MultiNode<typename T::PT, 3, typename T::ST, typename T::AT>, T> ||
    std::is_base_of_v<
        MultiNode<typename T::PT, 6, typename T::ST, typename T::AT>, T>;

template <typename T>
concept IsOrthTree = requires(T t) {
  { t.OrthTreeTag() } -> std::same_as<void>;
};

template <typename T>
concept IsKdTree = requires(T t) {
  { t.KdTreeTag() } -> std::same_as<void>;
};

template <typename T>
concept SupportsForceParallel = requires(T t) {
  { t.ForceParallel() } -> std::same_as<bool>;
};

template <typename T>
concept IsMaxStretchSplit = requires(T t) {
  { t.MaxStretchTag() } -> std::same_as<void>;
};

template <typename T>
concept IsRotateDimSplit = requires(T t) {
  { t.RotateDimTag() } -> std::same_as<void>;
};

struct FullCoveredTag {};
struct PartialCoverTag {};
}  // namespace cpdd
