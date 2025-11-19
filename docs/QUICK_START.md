# Quick Start

## With Docker
```bash
./docker-run.sh build
./docker-run.sh run --data-path /your/data/path
```

See [Docker manual](DOCKER_QUICK_REFERENCE.md) for complete details.

## Manual Setup
### Requirements

Necessary:

- `cmake` >= 3.15
- `g++` or `clang++` with C++20 features support (Tested with g++14 and clang-19) on Linux machines.

Optional:

- [jemalloc](https://github.com/jemalloc/jemalloc), slightly memory allocation improvement.
- [NUMA control](https://manpages.ubuntu.com/manpages/trusty/man8/numactl.8.html), improve the performance for parallelism.

## Getting Code
Try:

1. Clone the repository.

```bash
git clone git@github.com:ucrparlay/SpaceTreeLib.git
cd SpaceTreeLib
```

2. Initialize the submodule:

```bash
git submodule update --init
```

3. Compilation (with `Release` and `jemalloc` disabled)
```bash
mkdir build && cd build
cmake -DDEBUG=OFF -DJEMA=OFF ..
make -j
```

4. Run some toy examples:
```bash
./example/run_examples.sh
```

For more detailed usage of the PSI library, please checkout the [manual](MANUAL.md).
