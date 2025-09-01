# PSI: A library for Parallel Spatial Partition Trees

## Requirements

Necessary:

- CMake >= 3.15
- g++ or clang++ with C++20 features support (Tested with g++13 and clang17) on Linux machines, we suggest clang++ for better performance.
- We use [ParlayLib](https://github.com/cmuparlay/parlaylib) to support fork-join parallelism and some parallel primitives. It is provided as a submodule in our repository.

Optional:

- [jemalloc](https://github.com/jemalloc/jemalloc), slightly memory allocation improvement.
- [NUMA control](https://manpages.ubuntu.com/manpages/trusty/man8/numactl.8.html), improve the performance for parallelism.

## Getting Code

Try:

1. Clone the repository.

```{bash}
git clone git@github.com:ucrparlay/SpaceTreeLib.git
cd SpaceTreeLib
```

2. Initialize the submodule:

```bash
git submodule update --init
```

File structure:

```bash
├── benchmark
├── include
│   ├── parlaylib
│   └── psi
│       ├── base_tree_impl
│       │   └── tree_op
│       ├── cover_tree_impl
│       ├── dependence
│       │   ├── cpam
│       │   └── space_filling_curve
│       ├── kd_tree_impl
│       ├── orth_tree_impl
│       ├── p_tree_impl
│       └── r_tree_impl
├── script
├── static_analysis
└── tests
```

|    Name     |          Usage          |
| :---------: | :---------------------: |
| `benchmark` | Stores sample benchmark |
|  `include/psi`  |  Source of `PSI`   |
| `include/parlaylib` |   Provide parallelism   |
|  `scipts`   | Scripts for experiments |
|   `tests`   |   Helpers for testing   |

## Compilation

```bash
mkdir build && cd build
cmake -DDEBUG=OFF ..
make
```

Useful flags:

|   Name   |                 Usage                  | Default Value |
| :------: | :------------------------------------: | :-----------: |
| `DEBUG`  |         Compile in debug mode          |     `OFF`     |
| `SERIAL` |        Disable the parallelism         |     `OFF`     |
|  `CGAL`  | Compile all executable related to CGAL |     `OFF`     |
|  `JEMA`  |  Allocate the memory using `jemalloc`  |     `OFF`     |

More options can be found in `CMakeLists.txt`.

## Usage

Implemented in `tests/test.cpp`. See also `tests/cgal.cpp` and `tests/zdtree/neighborsTime.C`.

### Command line

```bash
./test -p [input_path] -d [dimension] -t [batch_mode] -r [rounds] -q [query_type] -i [read_insert_file_flag]
```

Parameters:

| Name |                       Usage                       |          Sample          |
| :--: | :-----------------------------------------------: | :----------------------: |
| `-p` |                 Input points path                 | `-p benchmark/sample.in` |
| `-d` |                Dimension of points                |          `-d 3`          |
| `-t` |  Batch mode, see [below](#test-framework-format)  |          `-t 0`          |
| `-r` |         How many rounds one test case run         |          `-r 3`          |
| `-q` | Query type, see [ below ](#test-framework-format) |          `-q 1`          |
| `-i` |     Whether to read the file for batch update     |          `-i 0`          |

### Example

Under `build/` folder, first generate two inputs, $P_1$ for tree construction and $P_2$ for batch update:

```{bash}
./data_generator ../benchmark 10000 3 2 0
```

which will generate two files named `1.in` and `2.in` in `benchmark/10000_3/`. Then to build the tree $P_1$, insert $P_2$, after which perform the $10$-NN search, try:

```{bash}
./test -p ../benchmark/10000_3/1.in -d 3 -t 1 -q 1 -r 3 -i 1
```

To parse the output, see [ Test Framework Format ](#test-framework-format) below.

### Default setting

In default, the PSIs stores all coordinates of points in 64-bit integer (`long`). The balancing parameter is set to $0.3$, the leaf wrap is $32$. It builds $6$ levels of tree at once. Different values on different machine may influence the performance dramatically. See our paper for more explanation.

## Test Framework Format

Implemented in `tests/test_framework.h`.

## Graph Generator
PSI provides two types of parallel data generator: `Uniform` and `Varden`.

Usage:
```{bash}
# In the folder build
make data_generator
./data_generator -p [output_path] -d [dimension] -n [points_num] -f [files_num] -t [0:uniform, 1:varden] 
```
It will create folders named as `uniform_bigint/` or `ss_varden_bigint/` under the directory specified by `output_path/`. The output file is numbered from `1.in` to `files_num`.

The data file begins with two integer, namely `points_num` and the `dimension`, follows by `points_num` number of lines, each line contains the coordinates for each point, separated by space.



