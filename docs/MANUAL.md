# PSI User Manual

PSI provides three parallel spatial partition trees:

| Header | Paper name | What it is |
|---|---|---|
| `psi/kd_tree.h` | Pkd-tree | parallel binary $k$d-tree |
| `psi/orth_tree.h` | P-Orth tree | parallel quad/oct-tree, any dimension |
| `psi/p_tree.h` | SPaC-tree | parallel 1-D tree over a Morton or Hilbert code |

There is also an AI-generated
[wiki](https://deepwiki.com/ucrparlay/SpaceTreeLib) covering how things work
internally. Contributions welcome 🤗.

Every code block below is compiled by `tests/unit/test_manual.cpp`, so if the
API moves and this file is not updated, the build breaks.

## Use PSI in your project

PSI is header-only. It needs C++20, parlaylib and libmorton.

```cmake
# Already have the repository, e.g. as a submodule:
add_subdirectory(SpaceTreeLib)
target_link_libraries(myapp PRIVATE PSI::PSI)
```

```cmake
# Or install it once (cmake --install) and find it:
find_package(PSI REQUIRED)
target_link_libraries(myapp PRIVATE PSI::PSI)
```

Directly, without CMake:

```bash
g++ -std=c++20 -O3 -pthread myapp.cpp \
    -I SpaceTreeLib/include \
    -I SpaceTreeLib/include/parlaylib/include \
    -I SpaceTreeLib/include/libmorton/include
```

`include/` is the only PSI include root, and every header is reached as
`psi/…`. Build your own targets with optimisation: PSI's `-O3` applies to
PSI's own targets, not to yours, and an unoptimised spatial index looks broken
rather than slow. Add `-march=native` if you are measuring.

Build options worth knowing:

| Option | Default | Meaning |
|---|---|---|
| `PSI_NATIVE_ARCH` | `ON` | `-march=native -mcx16`. Turn off for a binary that runs on another machine. |
| `PSI_ASSERTS` | `OFF` | keep `assert()` in an optimised build; the configuration to develop against |
| `PSI_BUILD_TESTS` | on if top-level | the `ctest` suite |
| `PSI_BUILD_BENCHMARKS` | on if top-level | the drivers and baselines; also what pulls in Boost |
| `PSI_BUILD_EXAMPLES` | on if top-level | the three example programs |
| `DEBUG` | `OFF` | unoptimised build with asserts and `validate()` |
| `CGAL` | `OFF` | adds the `kd_ccp` / `p_ccp` correctness oracles |

## Step-by-step guide

Taking the Pkd-tree; the other two differ only where noted under
[Choosing and configuring a tree](#choosing-and-configuring-a-tree). The
complete programs are in [example/](../example/).

### 1. A point type

```cpp
#include "psi/kd_tree.h"
#include "psi/dependence/splitter.h"

using coord_type = long;
using point = psi::basic_point<coord_type, 2>;
```

`coord_type` may be an integer or floating-point type. The dimension is a
template argument, so it is fixed at compile time.

### 2. A split rule

A split rule is a dimension rule paired with a partition rule:

```cpp
using split_rule = psi::orthogonal_split_rule<psi::rotate_dim<point>,
                                              psi::spatial_median<point>>;
```

| Dimension rule | Picks the split dimension by |
|---|---|
| `psi::rotate_dim<Point>` | cycling through the dimensions by depth |
| `psi::max_stretch_dim<Point>` | the widest side of the bounding box |

| Partition rule | Splits at |
|---|---|
| `psi::object_median<Point>` | the median point, so the two sides hold equal counts |
| `psi::spatial_median<Point>` | the middle of the box, so the two sides are equal in space |

### 3. A traits type

A traits type carries everything a tree is configured with: the point type,
the split rule, the augmentations, the skeleton height and the imbalance
ratio. Name it once and hand it to any tree that suits those points.

```cpp
using traits = psi::tree_traits<point, split_rule>;
```

Only the first two have no default. `psi/tree_traits.h` lists the rest.

### 4. Build

```cpp
psi::kd_tree<traits> tree;

parlay::sequence<point> points(1000);
/* ... fill points ... */
tree.build(parlay::make_slice(points));
```

`build` copies what it is given, so the input may go out of scope afterwards.
Calling `build` again on a live tree replaces its contents.

### 5. Query

Every point inside a box, and the `k` nearest neighbours of a point:

```cpp
typename decltype(tree)::box_type box;
box.first[0] = 400;  box.first[1] = 400;    /* lower corner */
box.second[0] = 600; box.second[1] = 600;   /* upper corner */

parlay::sequence<point> inside = tree.range_query(box);
auto nearest = tree.knn(query_point, 10);   /* sorted, nearest first */
```

`knn` gives back `(point, squared distance)` pairs by value, so they stay
valid after the tree changes. Both return fewer results than asked for when
the tree holds fewer points, and both return nothing on an empty tree.

That is all most callers need. The forms below let you own the storage
instead, which is what the benchmarks use:

```cpp
auto [count, count_log] = tree.range_count(box);

parlay::sequence<point> found(tree.get_size());
auto [written, query_log] = tree.range_query(box, parlay::make_slice(found));
/* found[0 .. written) are the points inside box */
```

`range_query` writes into the buffer you supply and returns how many points it
wrote. A buffer smaller than the answer is not an error: the result is
truncated and `written` tells you so. The second value in each pair is a
counter block for profiling; ignore it unless you are measuring.

```cpp
using dis_type = typename point::dis_type;
using nn_pair = std::pair<std::reference_wrapper<point>, dis_type>;

size_t const k = 10;
parlay::sequence<nn_pair> result(k, nn_pair(std::ref(points[0]), 0));
psi::bounded_queue<point, nn_pair> queue(parlay::make_slice(result));

tree.knn(query_point, queue);
/* result[0 .. queue.size()) hold the neighbours, unsorted, with squared
 * distances. queue.size() is less than k when the tree holds fewer points. */
```

The queue must be pre-filled, because `std::reference_wrapper` has no default
constructor; any point will do as the filler. Unlike `knn(q, k)`, these
results *reference* storage owned by the tree, so do not use them after
modifying it.

### 6. Update

```cpp
tree.batch_insert(parlay::make_slice(more_points));
tree.batch_delete(parlay::make_slice(points_to_remove));  /* must be present */
tree.batch_diff(parlay::make_slice(maybe_present));       /* tolerates absent */
```

`batch_delete` requires every point to be in the tree; `batch_diff` is the one
to use when that may not hold. The tree rebuilds imbalanced subtrees on its
own as updates accumulate.

### 7. Read everything back

```cpp
parlay::sequence<point> all(tree.get_size());
size_t n = tree.flatten(parlay::make_slice(all));
```

`flatten` needs a buffer of exactly `get_size()`; given any other size it
writes nothing and returns 0.

### 8. Clean up

The destructor frees the tree. `delete_tree()` does it early and may be
called more than once. A tree is move-only: it owns its nodes, so copying it
would free them twice.

## Choosing and configuring a tree

[Which tree do I want?](TREE_ANSWER.md) compares them. Concretely:

```cpp
/* Pkd-tree: any dimension rule, any partition rule */
psi::kd_tree<psi::tree_traits<point, split_rule>> kd;

/* P-Orth tree: partition rule must be spatial_median, and the skeleton
 * height must be a multiple of the dimension. Both are static_asserts;
 * orth_tree_traits picks a height that fits. The arity follows the point. */
psi::orth_tree<psi::orth_tree_traits<point, split_rule>> orth;

/* SPaC-tree: the point must carry the curve code, so it is an aug_point */
struct curve_code {
        using id_type = int_fast32_t;
        using curve_code_type = uint64_t;
        curve_code() : code(0), id(0) {}
        void set_member(curve_code_type const &v) { code = v; }
        bool operator<(curve_code const &r) const
        {
                return code == r.code ? id < r.id : code < r.code;
        }
        bool operator==(curve_code const &r) const { return id == r.id; }
        friend std::ostream &operator<<(std::ostream &o, curve_code const &r)
        {
                return o << r.code << " " << r.id;
        }
        curve_code_type code;
        id_type id;
};
using aug_pt = psi::aug_point<coord_type, 2, curve_code>;
using curve = psi::spatial_filling_curve<psi::morton_curve<aug_pt>>;
psi::p_tree<psi::tree_traits<aug_pt, curve>> p;
```

`psi::hilbert_curve<Point>` is the other curve. For `p_tree`, two points are
the same point when their `id`s match, so give every point a distinct one.

## Augmentations

Each node carries an augmentation. By default it is a bounding box, and you
need to write nothing:

```cpp
psi::kd_tree<psi::tree_traits<point, split_rule>> tree;
/* the traits default to box_leaf_aug and box_interior_aug */
```

Supply your own only if a bounding box is not what you want to store. The
contracts are the two concepts in `psi/dependence/concepts.h`:

```cpp
template <typename A, typename Slice>
concept leaf_augmentation = requires(A a, Slice in) {
        A(); A(in); a.update_aug(in); a.reset();
};

template <typename A>
concept interior_augmentation = requires(A a, bool flag, size_t n) {
        a.set_parallel_flag(flag);
        a.reset_parallel_flag();
        { a.get_parallel_flag_ini_status() } -> std::convertible_to<bool>;
        { a.force_parallel(n) } -> std::convertible_to<bool>;
};
```

Two things the concepts cannot state. First, an interior augmentation also
needs `create` and `update` member templates, because they take the very
`interior_type` that holds the augmentation, which is not nameable where the
concept is checked — copy their shape from `psi/dependence/augmentation.h`.
Second, the library looks for a **`get_box()`** member to decide whether a
node knows its own bounds; an augmentation that stores a box under a different
name silently falls back to the slower hyperplane paths.

`psi::box_leaf_aug` and `psi::box_interior_aug` in
`psi/dependence/augmentation.h` are the worked example.

## Parallelism

`build`, the batch updates and `flatten` are internally parallel. **A single
query is not**: `knn`, `range_count` and `range_query` each run on one thread,
so a loop over query points uses one core. Parallelise the batch yourself:

```cpp
parlay::parallel_for(0, queries.size(), [&](size_t i) {
        /* each iteration does its own query */
});
```

The thread count comes from parlaylib: set `PARLAY_NUM_THREADS`, or build with
`-DSERIAL=ON` for a single-threaded library. Queries do not mutate the tree, so
running many at once is safe; updates are not concurrent with anything.

## API reference

Common to all three trees. `Range` is any random-access range of points.

| Call | Returns | Notes |
|---|---|---|
| `build(Range&&)` | — | copies the input; replaces any existing contents |
| `batch_insert(Range&&)` | — | |
| `batch_delete(Range&&)` | — | every point must be present |
| `batch_diff(Range&&)` | — | tolerates points that are absent |
| `flatten(Range&& out)` | `size_t` | needs `out.size() == get_size()`, else writes nothing and returns 0 |
| `knn(Point const&, size_t k)` | `sequence<pair<Point, dis_type>>` | allocates; sorted, nearest first; points by value |
| `range_query(box_type const&)` | `points_type` | allocates and returns what it found |
| `knn(Point const&, bounded_queue&)` | profiling counters | results land in the queue; `queue.size()` is how many |
| `range_count(box_type const&)` | `pair<size_t, counters>` | |
| `range_query(box_type const&, Range&& out)` | `pair<size_t, counters>` | truncates to `out.size()` |
| `get_size()` | `size_t` | number of points |
| `empty()` | `bool` | |
| `delete_tree()` | — | idempotent |
| `get_tree_name()` | `char const*` | parsed by the plotting scripts; do not change |

`orth_tree::build` takes an optional bounding box, and `set_bounding_box` pins
one before building — useful when later inserts fall outside the initial
extent, since an orth tree does not grow its box. Inserting outside it is an
error.

Member types you will use: `point_type`, `box_type`, `points_type`,
`slice_type`, `coord_type`, `dis_type`. All of them come from the traits, so
`traits::box_type` and `tree::box_type` are the same type -- name whichever
reads better where you are.

## Extending PSI

Three extension points, in increasing order of effort.

### A new curve (p_tree)

The cheapest. Write a type with `encode()` in the shape of
`psi::morton_curve` / `psi::hilbert_curve` (`psi/dependence/splitter.h`) and
pass `psi::spatial_filling_curve<your_curve>` as `p_tree`'s split rule. The
code is applied during the build by `cpam_sample_sort`; `p_tree` itself never
calls the rule.

### A new split rule (kd_tree, orth_tree)

Derive from `psi::base_split_dim_rule` or `psi::base_split_partition_rule` and
implement every pure virtual — including `allow_rebuild()`, which batch update
asks before rebuilding a subtree.

Then add a branch for it in `orthogonal_split_rule::handle_undivided`, which
decides what happens when the chosen dimension cannot be split. Without one
you get a `static_assert` naming the rule; that assert is the whole reason it
is worth mentioning here.

### A new tree type

Not "derive and override". A derived tree has to supply, by name:

| What | Why |
|---|---|
| `leaf_type`, `interior_type` | every algorithm in `base_tree_impl/` takes them as template arguments |
| `splitter_type` | the build; the split rule itself comes from the traits |
| `node_arr_type`, `node_regions`, `splitter_num` | multi-way trees only; the skeleton builder reads them off the derived class |
| a tag method plus a matching `is_*` concept | how the shared code tells the families apart |
| `build_recursive` / `serial_build_recursive` | the split rule calls back into these, so the signature is fixed by the family |
| `friend split_rule_type` and the two `rebuild_*` friends | the callbacks reach private members |
| `delete_tree()` | pure virtual on the base |

`kd_tree` is the smaller model to copy; `orth_tree` shows the multi-way side.
Start by copying one wholesale and deleting, rather than deriving from scratch.

Two things the base cannot check for you: `validate()` has no branch for a node
type that is neither binary nor multi-way (it says so and continues), and the
augmentation concepts cannot state the `create`/`update` requirement.

## Checking a change

```bash
cmake -S . -B build -DPSI_ASSERTS=ON && cmake --build build -j
ctest --test-dir build              # the unit suite, well under a second
./example/run_examples.sh           # the examples check their own answers
```

For the CGAL oracle, configure with `-DCGAL=ON` and run
`script/checkCorrect.sh fast <data-root> build`. Generate a data root with the
tool below.

## A fast data generator

```bash
# In the build directory
./data_generator -p <out-dir> -n <points> -d <dim> -file_num 2 -varden 0 -seed 1
```

Files land in `<out-dir>/uniform/<n>_<d>/1.in` — or `uniform_bigint` when
`-axis_max` is above 10^6, which is the default, and `ss_varden…` with
`-varden 1`. Only dimensions 2, 3, 5, 7, 9, 12 and 16 are supported; anything else
throws. `-seed` makes the output reproducible; without it a fixed default
is used and printed.

## File organization

```
include/
├── baselines/          competitor implementations, benchmark only
│   ├── cpam_raw/       pristine CPAM and PAM
│   ├── zdtree/
│   └── zdtree_3d/
├── libmorton/          submodule, Morton codes
├── parlaylib/          submodule, parallel primitives
└── psi/
    ├── base_tree.h     shared machinery
    ├── base_tree_impl/
    ├── tree_traits.h   the configuration a tree takes, and the box arithmetic
    ├── tree_traits_impl/
    ├── dependence/     points, splitters, concepts, augmentations, cpam
    ├── kd_tree.h       + kd_tree_impl/
    ├── orth_tree.h     + orth_tree_impl/
    └── p_tree.h        + p_tree_impl/
example/                three worked programs, all of which check themselves
tests/unit/             the correctness suite
tests/                  benchmark drivers
script/, script_ae/     experiment drivers for the two papers
```

`psi/dependence/cpam/` is a fork of CPAM that `p_tree` is built on; see
[THIRD_PARTY.md](../THIRD_PARTY.md).
