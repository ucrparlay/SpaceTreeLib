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
- `g++` >= 14 or `clang++` >= 19, on Linux. Older versions are rejected at
  configure time.
- Boost headers (`libboost-dev`), for the benchmark drivers. Not needed if you
  only want the library: configure with `-DPSI_BUILD_BENCHMARKS=OFF`.

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

2. Initialize the submodules (the build stops with a clear message if you
   skip this):

```bash
git submodule update --init --recursive
```

3. Compilation (with `Release` and `jemalloc` disabled)
```bash
cmake -S . -B build -DDEBUG=OFF -DJEMA=OFF
cmake --build build -j
```

4. Check it works, from the repository root:
```bash
ctest --test-dir build      # the correctness suite, under a second
./example/run_examples.sh   # three programs that verify their own answers
```

For more detailed usage of the PSI library, please checkout the [manual](MANUAL.md).
