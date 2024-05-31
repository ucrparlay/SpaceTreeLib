#include <type_traits>
#include <tuple>

namespace cpdd {

// NOTE: the legacy iterator which only increase counter
template<typename T>
class counter_iterator {
   private:
    struct accept_any {
        template<typename U>
        accept_any& operator=(const U&) {
            return *this;
        }
    };

   public:
    typedef std::output_iterator_tag iterator_category;

    counter_iterator(T& counter) : counter(counter) {}
    counter_iterator(const counter_iterator& other) : counter(other.counter) {}

    bool operator==(const counter_iterator& rhs) const { return counter == rhs.counter; }
    bool operator!=(const counter_iterator& rhs) const { return counter != rhs.counter; }

    accept_any operator*() const {
        ++counter.get();
        return {};
    }

    counter_iterator& operator++() {  // ++a
        return *this;
    }
    counter_iterator operator++(int) {  // a++
        return *this;
    }

   protected:
    std::reference_wrapper<T> counter;
};

// NOTE: check whether the type is a pair
template<typename>
struct is_pair : std::false_type {};

template<typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {};

// NOTE: tag for different leaf allocations
struct AllocNormalLeafTag {};
struct AllocDummyLeafTag {};
struct alloc_fat_leaf_tag {};
struct alloc_thin_leaf_tag {};

}  // namespace cpdd
