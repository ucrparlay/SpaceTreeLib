# PSPT: A library for Parallel Spatial Partition Trees

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
.
├── benchmark
├── include
├── parlaylib
├── script
└── tests
```

|    Name     |          Usage          |
| :---------: | :---------------------: |
| `benchmark` | Stores sample benchmark |
|  `include`  |  Source of `PSPT`   |
| `parlaylib` |   Provide parallelism   |
|  `scipts`   | Scripts for experiments |
|   `tests`   |   Helpers for testing   |

## Compilation

```bash
mkdir build && cd build
cmake -DDEBUG=OFF ..
make
```

For better performance, please use `clang++` for compilation, i.e.,

```{bash}
cmake -DDEBUG=OFF -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ..
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

In default, the PSPTs stores all coordinates of points in 64-bit integer (`long`). The balancing parameter is set to $0.3$, the leaf wrap is $32$. It builds $6$ levels of tree at once. Different values on different machine may influence the performance dramatically. See our paper for more explanation.

## Test Framework Format

Implemented in `tests/test_framework.h`.




