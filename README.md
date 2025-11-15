# PSI: A library for Parallel Spatial Indexes

## Docker Setup (Recommended for Artifact Evaluation)

**Hardware Requirements:**
- **Memory**: 512GB RAM for full experiments
- **CPU**: All available cores
- **Disk**: 100GB+ free space

**Quick Start:**
```bash
# Build the Docker image
./docker-run.sh build

# Run experiments (specify where to store generated data)
./docker-run.sh run --data-path /mnt/large-disk/data

# Or run full evaluation
./docker-run.sh full --data-path /mnt/large-disk/data
```

**Parameters:**
- `--data-path PATH` - Host directory for generated datasets (required)
- `--node-size SIZE` - Number of data points (default: 1000000000)
- `--memory SIZE` - Memory limit (e.g., 512g)
- `--cpus NUM` - CPU limit (default: all cores)

**Results:** Found in `results/`, `logs/`, `plots/` directories on host.

**Documentation:** See [doc/](doc/) folder for complete guides.

## Requirements

Necessary:

- CMake >= 3.15
- g++ or clang++ with C++20 features support (Tested with g++14 and clang-19) on Linux machines.
- We use [ParlayLib](https://github.com/cmuparlay/parlaylib) to support fork-join parallelism and some parallel primitives. It is provided as a submodule in our repository.

Optional:

- [jemalloc](https://github.com/jemalloc/jemalloc), slightly memory allocation improvement.
- [NUMA control](https://manpages.ubuntu.com/manpages/trusty/man8/numactl.8.html), improve the performance for parallelism.

## Getting Code

Try:

1. Clone the repository.

```{bash}
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
| `include/baselines` |   implmentation of baselines   |
|  `scipts`   | Scripts for experiments |
|   `tests`   |   Helpers for testing   |

## Compilation

```bash
mkdir build && cd build
cmake -DDEBUG=OFF -DCGAL=OFF -DJEMA=OFF ..
make
```

Useful flags:

|   Name   |                 Usage                  | Default Value |
| :------: | :------------------------------------: | :-----------: |
| `DEBUG`  |         Compile in debug mode          |     `OFF`     |
| `SERIAL` |        Disable the parallelism         |     `OFF`     |
|  `CGAL`  | used to check the correctness |     `OFF`     |
|  `JEMA`  |  Allocate the memory using `jemalloc`  |     `OFF`     |

More options can be found in `CMakeLists.txt`.

## Usage

Implemented in `tests/*_test.cpp`. 

### Command line

```bash
./${solver}$ -p [input_path] -d [dimension] -t [batch_mode] -r [rounds] -q [query_type] -i [read_insert_file_flag]
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
| `-T` |     The tree type     |          `-T 0`          |
| `-l` |     The split for the tree     |          `-l 0`          |

### Example

Under `build/` folder, first generate two inputs, $P_1$ for tree construction and $P_2$ for batch update:

```{bash}
./data_generator -p ../benchmark -n 10000 -d 2 -file_num 2 -varden 0
```

which will generate two files named `1.in` and `2.in` in `benchmark/10000_3/`. Then to build a quad-tree $P_1$, insert $P_2$, after which perform the $10$-NN search, try:

```{bash}
./kd_test -p ../benchmark/uniform/10000_2/1.in -d 2 -t 1 -q 1 -r 3 -i 1 -T 1 -l 3 
```

To parse the output, see [ Test Framework Format ](#test-framework-format) below.

## Test Framework Format

Implemented in `tests/test_framework.h`.

## Graph Generator
PSI provides two types of parallel data generator: `Uniform` and `Varden`.

Usage:
```{bash}
# In the folder build
make data_generator
./data_generator -p [output_path] -d [dimension] -n [points_num] -file_num [files_num] -varden [0:uniform, 1:varden] 
```
It will create folders named as `uniform_bigint/` or `ss_varden_bigint/` under the directory specified by `output_path/`. The output file is numbered from `1.in` to `files_num`.

The data file begins with two integer, namely `points_num` and the `dimension`, follows by `points_num` number of lines, each line contains the coordinates for each point, separated by space.



