# PSI User Manual
The PSI library provides three parallel spatial partition trees that are widely used, namely:
- Parallel kd-tree (Pkd-tree);
- Parallel quad/oct-tree (POrth-tree);
- Parallel 1d-tree based on spatial filling curves(SPaC-tree);

Thanks for AI, there is a generally good [wiki](https://deepwiki.com/ucrparlay/SpaceTreeLib) for this library, where you can find out how things are operated, and maybe get some start ideas if you want to contribute some code ☺️.

Any contribution is welcomed 🤗!

## File Organization
```{bash}
include/
├── baselines
│   ├── boost_rtree
│   ├── cpam_raw
│   ├── zdtree
│   └── zdtree_3d
├── libmorton
├── parlaylib
└── psi
    ├── base_tree_impl
    ├── base_tree.h
    ├── dependence
    ├── kd_tree_impl
    ├── kd_tree.h
    ├── orth_tree_impl
    ├── orth_tree.h
    ├── p_tree_impl
    ├── p_tree.h
```

- `baselines/`: the codes for baselines.
- `libmorton/`: the [open-source library](https://github.com/Forceflow/libmorton) to compute the Morton/Z code quickly.
- `parlaylib/`: a [toolkit for parallel algorithms](https://github.com/cmuparlay/parlaylib).
- `psi/`: the main sources for the PSI library.
- `psi/base_tree.h`: the underlying tree structure provides type definition, geometry toolkit and general APIs.
- `psi/kd_tree.h`: known as the **Pkd-tree** the implementation for the parallel $k$d-tree.
- `psi/orth_tree.h`: known as the **P-Orth tree**, the implementation for the parallel Quad/Orth-tree (can be extened to arbitrary higher dimensions).
- `psi/p_tree.h`: known as the **SPaC-tree**, the implementation for the parallel 1-D spatial trees based on the Morton/Z code or Hilbert code.

## Step-by-step guide
We will take the **Pkd-tree** as an example, other trees are similar.
All source files for examples can be found in the [example/](../example/).

Everything start with the definition of the `Point`:
```{c++}
// Type for each coordinate
using Coord = long;

// Define augmentation structure for points (stores an ID)
// WARN: All functions must be defined
struct AugId {
  using IdType = int;
  IdType id;

  bool operator<(AugId const& rhs) const { return id < rhs.id; }
  bool operator==(AugId const& rhs) const { return id == rhs.id; }
  friend std::ostream& operator<<(std::ostream& os, AugId const& rhs) {
    os << rhs.id;
    return os;
  }
};

// Define point type: 2D points with augmented ID
using Point = psi::AugPoint<Coord, 2, AugId>;
```

With `Point`, we can now use all functionalities provided in `Basetree`, e.g.,
```{c++}
using BT = psi::BaseTree<Point>;
auto box = BT::GetBox(input); // demo code, to get the bounding box for the input
```

You want to augment some info on the tree nodes, both the leaf nodes, and the interior nodes, you can define the augmentation in this way:
<details>
<summary>click to expand</summary>

```{c++}
// Leaf augmentation: stores bounding box
// WARN: All functions must be defined
template <class BaseTree>
struct LeafAugBox {
  using Box = typename BaseTree::Box;      // The bounding box type
  using Slice = typename BaseTree::Slice;  // parallel data container, similar
                                           // to std::ranges

  // Constructors
  LeafAugBox() : box(BaseTree::GetEmptyBox()) {}
  LeafAugBox(Box const& _box) : box(_box) {}
  LeafAugBox(Slice In) : box(BaseTree::GetBox(In)) {}

  // Return the bounding box
  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  // Update the augmentation information
  void UpdateAug(Slice In) { this->box = BaseTree::GetBox(In); }

  // Reset the augmentation information
  void Reset() { this->box = BaseTree::GetEmptyBox(); }

  Box box;
};

// Interior node augmentation: stores bounding box and parallel build flag
// WARN: All functions must be defined
template <class BaseTree>
struct InteriorAugBox {
  using Box = typename BaseTree::Box;

  // Constructors
  InteriorAugBox() : box(BaseTree::GetEmptyBox()) {
    force_par_indicator.reset();
  }
  InteriorAugBox(Box const& _box) : box(_box) { force_par_indicator.reset(); }

  // Get the bounding box
  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  // Given two child nodes, create the augmentation value (bounding box) for
  // the interior node
  template <typename Leaf, typename Interior>
  static Box Create(psi::Node* l, psi::Node* r) {
    return BaseTree::GetBox(BaseTree::template RetriveBox<Leaf, Interior>(l),
                            BaseTree::template RetriveBox<Leaf, Interior>(r));
  }

  // Update the augmentation information for the interior node
  template <typename Leaf, typename Interior>
  void Update(psi::Node* l, psi::Node* r) {
    this->box = Create<Leaf, Interior>(l, r);
  }

  // Below are required to ensure the granularity of parallelism
  // Mark this node to be update in parallel or not in following operations
  void SetParallelFlag(bool const flag) {
    this->force_par_indicator.emplace(flag);
  }

  // Reset the parallel flag
  void ResetParallelFlag() { this->force_par_indicator.reset(); }

  // Get the initial status of the parallel flag
  bool GetParallelFlagIniStatus() {
    return this->force_par_indicator.has_value();
  }

  // Whether to force parallel update for this node
  bool ForceParallel(size_t sz) const {
    return this->force_par_indicator.has_value()
               ? this->force_par_indicator.value()
               : sz > BaseTree::kSerialBuildCutoff;
  }

  // Reset the augmentation information and parallel flag
  void Reset() {
    this->force_par_indicator.reset();
    this->box = BaseTree::GetEmptyBox();
  }

  Box box;
  std::optional<bool> force_par_indicator;
};
```
</details>

Of course, as a spatial partition tree, we can choose how to split the space:
```{c++}
// Define split rule: max stretch dimension + object median
using SplitRule = psi::OrthogonalSplitRule<psi::MaxStretchDim<Point>,
                                           psi::ObjectMedian<Point>>;

// Alternative split rule: rotate dimension + spatial median
using AnotherSplitRule =
    psi::OrthogonalSplitRule<psi::RotateDim<Point>, psi::SpatialMedian<Point>>;
```

Now we can define the tree type with all building blocks above:
```{c++}
// Define KdTree type
using Tree = psi::KdTree<Point, SplitRule, LeafAugBox<BT>, InteriorAugBox<BT>>;
Tree tree;
```

Then we can use the tree as follows:
- Build the tree:
```{c++}
Points points;
auto points_copy = points;
tree.Build(points_copy);
```
- Batch Insert:
```{c++}
Points insert_points;
auto insert_copy = insert_points;
tree.BatchInsert(insert_copy);
```

- Batch Delete (assumes all points to be deleted are in the tree, use `BatchDiff` if you are not sure):
```{c++}
Points delete_points;
auto delete_copy = delete_points;
tree.BatchDelete(delete_copy);
// tree.BatchDiff(delete_copy);
```

- KNN query
```{c++}
int K = 10; // K for KNN
Point query_point; // Points to be queried

using DisType = typename Point::DisType; // the distance type we use
using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>; // how KNN candidates are represented in the output

parlay::sequence<nn_pair> knn_result(K, nn_pair(std::ref(points[0]), 0));// the output array
psi::kBoundedQueue<Point, nn_pair> bq(parlay::make_slice(knn_result));// a fast heap for query

auto* root = tree.GetRoot();
tree.KNN(root, query_point, bq); // do the query
```

- Range count and Range query
```{c++}
typename Tree::Box query_box; // a pair of point
Points range_result(n);  // Allocate max possible size
auto [count, logger] =
      tree.RangeQuery(query_box, parlay::make_slice(range_result));
```

A comprehensive example can be found [here](../example/kd_tree.h).

## A fast data generator
PSI also shipped a parallel data_generator for two distribution of points, namely, the `Uniform` (uniformly spread across a cube) and `Varden` (very skewed).

Usage:
```{bash}
# In the build
make data_generator
./data_generator -p [output_path] -d [dimension] -n [points_num] -file_num [files_num] -varden [0:uniform, 1:varden] 
```

It will create folders named as `uniform_bigint/` or `ss_varden_bigint/` under the directory specified by `output_path/`. The output file is numbered from `1.in` to `files_num`.

The data file begins with two integer, namely `points_num` and the `dimension`, follows by `points_num` number of lines, each line contains the coordinates for each point, separated by space.

