#ifndef PSPT_DEPENDENCE_CONCEPTS_H_
#define PSPT_DEPENDENCE_CONCEPTS_H_

#include <concepts>
#include <type_traits>
#include <utility>
#include <variant>

namespace pspt {

// NOTE: In memory of SFINE
// template<typename>
// struct is_pair : std::false_type {};
// template<typename T, typename U>
// struct is_pair<std::pair<T, U>> : std::true_type {};
// template<typename T>
// concept IsPair = is_pair<T>::value;

template <typename T>
concept IsPointer = std::is_pointer_v<T>;

template <typename T, typename Node>
concept CheckPointerToNode =
    std::is_pointer_v<T> && std::is_base_of_v<Node, std::remove_pointer_t<T>>;

template <typename T>
concept IsPair = requires {
  requires std::same_as<
      T, std::pair<typename T::first_type, typename T::second_type>>;
};

template <typename T, typename Point>
concept IsBox = requires {
  requires IsPair<T> &&
               std::same_as<typename T::first_type, typename Point::BP> &&
               std::same_as<typename T::second_type, typename Point::BP>;
};

// NOTE: tag whether a node has non-trivial augmentation
template <typename T>
concept NodeHasNonTrivialAug =
    (!std::is_empty_v<T>) && (!std::same_as<T, std::monostate>);

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
// if func accepts a parameter arg, then call it; otherwise return func()
// directly
template <typename ReturnType, typename Func, typename Arg>
ReturnType InvokeWithOptionalArg(Func&& func, Arg&& arg) {
  if constexpr (CallableWithArg<ReturnType, Func, Arg>) {
    return func(std::forward<Arg>(arg));
  } else {
    return func();
  }
}

// NOTE: define the what is a binary node
template <typename T,
          template <typename, typename, typename> typename BinaryNodeT>
concept CheckBinaryNode = std::is_base_of_v<
    BinaryNodeT<typename T::PT, typename T::ST, typename T::AT>, T>;

// NOTE: helper for decide a multi-way node
template <
    typename T,
    template <typename, uint_fast8_t, typename, typename> typename MultiNodeT,
    uint_fast8_t... Ns>
concept IsMultiNodeHelper =
    (std::is_base_of_v<
         MultiNodeT<typename T::PT, Ns, typename T::ST, typename T::AT>, T> ||
     ...);

// NOTE: define the what is a multi-way node
template <typename T, template <typename, uint_fast8_t, typename,
                                typename> typename MultiNodeT>
concept CheckMultiNode =
    IsMultiNodeHelper<T, MultiNodeT, 2, 3, 4, 5, 6, 7, 8, 9, 10>;

// NOTE: define the what is a dynamic node
template <typename T,
          template <typename, typename, typename> typename DynamicNodeT>
concept CheckDynamicNode = std::is_base_of_v<
    DynamicNodeT<typename T::PT, typename T::ST, typename T::AT>, T>;

// NOTE: tag aug point
template <typename T>
concept IsAugPoint = requires(T t) {
  { t.AugPointTag() } -> std::same_as<void>;
};

// NOTE: tag contain bounding box
template <typename T>
concept HasBox = requires(T t) {
  { t.box };
};

// NOTE: tag a orth tree
template <typename T>
concept IsOrthTree = requires(T t) {
  { t.OrthTreeTag() } -> std::same_as<void>;
};

// NOTE: tag a kd tree
template <typename T>
concept IsKdTree = requires(T t) {
  { t.KdTreeTag() } -> std::same_as<void>;
};

// NOTE: tag a rtree
template <typename T>
concept IsRTree = requires(T t) {
  { t.RTreeTag() } -> std::same_as<void>;
};

// NOTE: tag a rtree
template <typename T>
concept IsCoverTree = requires(T t) {
  { t.CoverTreeTag() } -> std::same_as<void>;
};

template <typename T>
concept IsPTree = requires(T t) {
  { t.PTreeTag() } -> std::same_as<void>;
};

template <typename T>
concept SupportsForceParallel = requires(T t) {
  { t.ForceParallel() } -> std::same_as<bool>;
};

template <typename T>
concept IsMaxStretchDim = requires(T t) {
  { t.MaxStretchTag() } -> std::same_as<void>;
};

template <typename T>
concept IsRotateDimSplit = requires(T t) {
  { t.RotateDimTag() } -> std::same_as<void>;
};

template <typename T>
concept IsObjectMedianSplit = requires(T t) {
  { t.ObjectMedianTag() } -> std::same_as<void>;
};

template <typename T>
concept IsSpatialMedianSplit = requires(T t) {
  { t.SpatialMedianTag() } -> std::same_as<void>;
};

template <typename T>
concept IsHilbertCode = requires(T t) {
  { t.HilbertTag() } -> std::same_as<void>;
};

template <typename T>
concept IsMortonCode = requires(T t) {
  { t.MortonTag() } -> std::same_as<void>;
};
struct FullCoveredTag {};
struct PartialCoverTag {};

}  // namespace pspt

#endif  // PSPT_DEPENDENCE_CONCEPTS_H_
