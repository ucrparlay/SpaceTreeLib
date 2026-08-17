# Contributing to PSI

Contributions are welcome. PSI backs two papers, so a few of its rules are
stricter than they would otherwise be; they are all here.

## Build and check

```bash
cmake -S . -B build -DPSI_ASSERTS=ON && cmake --build build -j
ctest --test-dir build --output-on-failure   # under a second
./example/run_examples.sh                    # the examples check themselves
```

`-DPSI_ASSERTS=ON` keeps `assert()` in an optimised build. That is the
configuration to develop against: `-DDEBUG=ON` also enables the asserts but is
unoptimised, so it is 10-50x slower.

For a stronger check, compare against CGAL:

```bash
cmake -S . -B build -DCGAL=ON -DPSI_ASSERTS=ON && cmake --build build -j
build/data_generator -p /tmp/psi-data/uniform -n 20000 -d 2 -file_num 2 -varden 0
cd script && PSI_NODES=20000 ./checkCorrect.sh fast /tmp/psi-data ../build
```

Anything that could move a measured number needs re-measuring, and the pull
request should say so. That includes changes that look like pure cleanup.

## Formatting

```bash
cmake --build build --target format
```

Whatever `.clang-format` produces: tabs, 8 wide, 80 columns, Linux braces.
The tree is formatted with **clang-format 19.1.0**, which CI installs with
`pip install clang-format==19.1.0` -- versions disagree about this code, so
use that one if `make format` produces churn.
Every file must be a clang-format fixpoint; CI checks it. `make format` covers
`include/psi`, `include/baselines`, `tests` and `example`, and must never
reach `include/parlaylib` or `include/libmorton`, which are submodules.

## Naming

As the standard library names things: `snake_case` for classes, functions,
variables and concepts; member type aliases ending in `_type`; `CamelCase`
for template parameters; trailing underscore on private data members;
`UPPER_CASE` for macros. Constants are `snake_case` with no prefix.

A few names stay `CamelCase` because someone else owns them — `Tree` is
CPAM's and CGAL's nested type, and `tests/kd_ccp.cpp` typedefs CGAL names onto
themselves.

## Configuration

A tree takes one template parameter, a traits type, and `base_tree` inherits
it:

```cpp
using traits = psi::tree_traits<point, split_rule>;
psi::kd_tree<traits> tree;
```

A new alias or constant goes in `psi/tree_traits.h` -- in `point_traits` if it
follows from the point type alone, in `tree_traits` otherwise. Then re-export
it from the block at the top of `base_tree.h`: the traits are a dependent
base, so a name missing from that block is invisible to the impl files.

## Comments

Few and short, `/* */`, and about *why*. If code needs a comment to say what
it does, rename something instead. No `NOTE:`/`WARN:`/`PERF:` tags, no banner
lines, no commented-out code — git remembers it.

Write the simple version. No abstraction until there are three callers, no
helper used once, no option nobody asked for.

## What not to change

The parallel algorithms are frozen. Their control flow, work and depth are
out of scope unless a change is explicitly asked for:

- the `*_recursive` functions in `kd_tree_impl/` and `orth_tree_impl/`
- the skeleton tag values `bucket_num+1/+2/+3` in `inner_tree.hpp`
- the force-parallel flag and the `serial_build_cutoff` fallback
- the granularity argument on every `parallel_for`
- `-march=native`, `-O3`, `-mcx16`

Also out of scope:

- `include/baselines/` exists to compare against. Do not refactor it; change
  it only where it must follow a rename to keep compiling.
- `include/psi/dependence/cpam/` is a permanent PSI fork of CPAM (upstream is
  no longer maintained). It is deferred, not untouchable: leave it alone
  unless asked. See [THIRD_PARTY.md](THIRD_PARTY.md).

Two consequences of that last one, so they are not rediscovered: the seven
`-Wmaybe-uninitialized` warnings in a release build all come from the fork and
are benign, and `p_tree` has no `validate()` because its node type is CPAM's —
it is covered by the black-box tests under `tests/unit` instead.

## Load-bearing strings

`get_tree_name()` and `check_has_box()` return strings that
`script_ae/merge_*.py` parses positionally to label the figures. Changing them
silently mislabels results. The same goes for the column order in
`print_tree_param()`.

## Errors

`assert()` is for internal invariants. A check on something the caller
supplied belongs at a public entry point, as a documented result the caller
can see in a release build — an empty query box returns 0, a too-small output
buffer truncates and reports what was written. Do not add an assert that a
release build silently drops when the condition is a caller's mistake.

## Tests

New behaviour needs a test under `tests/unit/`, which includes only `psi/…`
headers — never `tests/test_framework.h`, which costs about 70 s and 1.7 GB
per translation unit and pulls in the baselines. Add the file to
`tests/unit/CMakeLists.txt`; it is registered with CTest automatically.

Code in `docs/MANUAL.md` is compiled by `tests/unit/test_manual.cpp`. If you
change the API, change both.
