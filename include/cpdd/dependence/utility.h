#include <type_traits>
#include <tuple>

namespace cpdd {

// NOTE: check whether the type is a pair
template<typename>
struct is_pair : std::false_type {};

template<typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {};

// NOTE: tag for different leaf allocations
struct AllocNormalLeafTag {};
struct AllocDummyLeafTag {};
struct AllocEmptyLeafTag {};
struct BinaryInteriorTag {};
struct MultiWayInteriorTag {};
struct alloc_fat_leaf_tag {};
struct alloc_thin_leaf_tag {};

}  // namespace cpdd
