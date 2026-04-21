# SIMD Sample Sort Branch Notes

This document summarizes the changes on branch `sean/sort`, how to build the project with SIMD sample sort enabled, and how to run the current build-only benchmark path.

## What Changed

The core change is that the SPaC-tree / `PTree` build path no longer sorts CPAM entries in the old `(AugType, Point*)` form. It now materializes a compact key-value representation:

- key: `CurveCode` (`uint64_t`)
- value: point address stored as `uint64_t`

The main type switch happens in [include/psi/p_tree.h](/home/xwang605/SpaceTreeLib/include/psi/p_tree.h:66):

```cpp
using key_entry_pointer =
    parlay::KVPair<typename Point::AT::CurveCode, uint64_t>;
```

This enables the CPAM sort path to work directly on integer keys and integer payloads, which matches the SIMD sample sort interface much better than the previous pointer-pair form.

## CPAM Sort Path

The sort adapter in [include/psi/dependence/cpam/cpam_sample_sort.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/cpam_sample_sort.h:33) now does three things:

1. Detects whether SIMD sample sort is available through `PSI_USE_SIMD_SAMPLE_SORT`.
2. Materializes `(curve_code, address)` pairs before sorting.
3. Falls back to the original CPAM sample sort path if the SIMD path is not enabled.

Key additions:

- `make_cpam_output_entry(...)` generalizes output generation for plain entries, augmentation keys, pointer pairs, and integral key-value pairs.
- `cpam_sample_sort_simd_pair(...)` precomputes filling-curve codes and sorts the materialized array with the SIMD implementation.
- `cpam_sample_sort_materialized_pair(...)` provides a non-SIMD precomputed fallback for the same compact representation.

In short, the branch changes the cost structure from "compute/compare richer objects during sort" to "precompute sortable integer keys once, then sort compact records".

More precisely, the two routes are:

- `SIMD ON`
  `fill code + materialize KV -> SIMD sort -> build`
- `SIMD OFF`
  the legacy CPAM sample-sort path, where code filling is interleaved with the sort rather than separated as a clean standalone phase

So if you want a step model:

- `SIMD ON`: `fill -> sort -> build`
- `SIMD OFF`: closer to `fill+sort(interleaved) -> build`

## CPAM Build Adaptation

Because the payload is now stored as an integerized address, CPAM build helpers were updated to reconstruct entries through a shared accessor:

- [include/psi/dependence/cpam/basic_node_helpers.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/basic_node_helpers.h:11)
- [include/psi/dependence/cpam/sequence_ops.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/sequence_ops.h:721)
- [include/psi/dependence/cpam/build.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/build.h:42)

That keeps the rest of the tree build logic mostly unchanged while allowing both legacy pointer payloads and the new integral payload form.

## Test Framework Changes

The benchmark harness in [tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:73) was also adjusted to fit the new sort/build path:

- `AugCode` now stores only the filling-curve code, not a point id.
- Input and insert sets for `PTree` are deduplicated by coordinates before building.
- Verification helpers were changed so checks do not rely on `aug.id` being present.
- A new build-only mode was added, and it now uses a dedicated helper for pure
  build benchmarking:

```cpp
if (kTag == 0 && kQueryType == 0) {
  BenchmarkBuildOnly<Point, Tree>(wp, tree, kRounds);
  return;
}
```

This means `-t 0 -q 0` now runs tree construction only.

The important detail is that `BenchmarkBuildOnly(...)` no longer reuses the
generic `BuildTree(...)` flow. `BuildTree(...)` was designed for downstream
insert/query tests and therefore performed an extra non-measured build after
timing. For the build-only path we now run:

1. one warmup `build + delete`
2. three measured `build + delete` iterations
3. report the average of the measured iterations only

So the build-only path now matches the intended "drop the first run and average
the next three" workflow.

## Build Command

After submodules are available, a working SIMD build is:

```bash
cmake -S . -B build-simd -DUSE_SIMD_SAMPLE_SORT=ON
cmake --build build-simd -j
```

If you also want the CPAM build-stage timing prints, enable:

```bash
cmake -S . -B build-simd \
  -DUSE_SIMD_SAMPLE_SORT=ON \
  -DPRINT_CPAM_BUILD_TIMING=ON
cmake --build build-simd -j --target p_test data_generator
```

If you only want the relevant binaries:

```bash
cmake --build build-simd -j --target p_test data_generator
```

If you want jemalloc enabled at build time, the project option is:

```bash
cmake -S . -B build-simd \
  -DUSE_SIMD_SAMPLE_SORT=ON \
  -DJEMA=ON
```

So the malloc-related option you remembered is `JEMA`, which enables jemalloc in this project.

## Uniform 1e9 Data Generation

The generator is documented in [docs/MANUAL.md](/home/xwang605/SpaceTreeLib/docs/MANUAL.md:232) and implemented in [tests/data_generate.cpp](/home/xwang605/SpaceTreeLib/tests/data_generate.cpp:1).

To generate one uniform dataset with 1e9 points and no washing/sorting:

```bash
./build-simd/data_generator \
  -p /path/to/data \
  -n 1000000000 \
  -d 2 \
  -file_num 1 \
  -varden 0
```

This creates:

```text
/path/to/data/uniform_bigint/1000000000_2/1.in
```

For 3D with the smaller axis range used by the AE scripts:

```bash
./build-simd/data_generator \
  -p /path/to/data \
  -n 1000000000 \
  -d 3 \
  -file_num 1 \
  -varden 0 \
  -axis_max 1000000
```

## Build-Only Run Command

Your current framework change makes the SPaC-tree build-only run:

```bash
PARLAY_NUM_THREADS=192 numactl -i all ./build-simd/p_test \
  -p /path/to/data/uniform_bigint/1000000000_2/1.in \
  -d 2 \
  -r 1 \
  -t 0 \
  -q 0 \
  -i 0 \
  -T 2 \
  -l 2
```

Parameter notes:

- `-T 2`: choose `PTree`
- `-l 2`: Morton curve
- `-l 1`: Hilbert curve
- `-t 0 -q 0`: build only
- `-i 0`: do not read an insert file

## Timing Output

With `-DPRINT_CPAM_BUILD_TIMING=ON`, building a `PTree` from an empty map prints an extra line such as:

```text
[CPAM build] route=simd-precompute-fill-sort n=... fill=... sort=... build=... total=...
```

or:

```text
[CPAM build] route=legacy-interleaved-fill-sort n=... fill=interleaved-with-sort sort+fill=... build=... total=...
```

Meaning:

- `route` identifies the path used
- `fill` exists as a separate phase only for the SIMD / precompute route
- `sort` is pure sorting time
- `sort+fill` is used for the legacy route because those two are intertwined
- `build` is the tree-construction phase after sorting

## Current 1e9 Morton Results

Command used:

```bash
numactl -i all ./build-simd/p_test \
  -p /tmp/psi-data/uniform_bigint/1000000000_2/1.in \
  -d 2 \
  -r 1 \
  -t 0 \
  -q 0 \
  -i 0 \
  -T 2 \
  -l 2
```

Interpretation:

- the first `[CPAM build]` line is the warmup run
- the last three `[CPAM build]` lines are the measured runs
- the final scalar printed by `p_test` is the average of the measured runs in
  the build-only harness

Measured averages over the last three runs:

| Route | Fill / sort stage | Build stage | CPAM total |
| --- | ---: | ---: | ---: |
| Legacy (`legacy-interleaved-fill-sort`) | `1.811655s` (`sort+fill`) | `0.692736s` | `2.504391s` |
| SIMD (`simd-precompute-fill-sort`) | `0.356934s fill + 1.145278s sort = 1.502212s` | `0.685634s` | `2.187846s` |

Observed differences:

- `fill+sort` improved by `0.309443s`, from `1.811655s` to `1.502212s`
  (`1.206x`, `17.08%` faster)
- `build` improved by only `0.007102s`, from `0.692736s` to `0.685634s`
  (`1.03%` faster, effectively flat)
- CPAM `total` improved by `0.316546s`, from `2.504391s` to `2.187846s`
  (`1.145x`, `12.64%` faster)

So the SIMD branch is not really accelerating the tree-construction kernel
itself. The gain is concentrated in the pre-build preparation path: code fill,
materialization, and sorting.

For the final scalar printed by the harness:

- legacy build-only average: `2.85636`
- SIMD build-only average: `2.63982`

That is an end-to-end build-only improvement of `0.21654s`
(`1.082x`, `7.58%` faster). This number is smaller than the CPAM-stage speedup
because it also includes framework overhead outside the `[CPAM build]` trace.

## Branch Summary

At this point the branch changes fall into five buckets:

1. CPAM sort/data-path changes
   - materialized `(CurveCode, address)` representation for `PTree`
   - SIMD sample sort adapter and a scalar materialized fallback
   - route/timing tracing for empty-map builds
2. CPAM build compatibility changes
   - helper updates so both pointer payloads and integerized-address payloads
     can be reconstructed during build
3. Test harness changes
   - `AugCode`
   - coordinate dedup for `PTree`
   - verification helpers that no longer require `aug.id`
   - dedicated build-only benchmark flow with explicit warmup semantics
4. Build/config integration
   - `USE_SIMD_SAMPLE_SORT`
   - `PRINT_CPAM_BUILD_TIMING`
   - Highway integration under `include/SIMD-Sample-Sort`
   - `.gitignore` updated to ignore `build-simd/`
5. Documentation
   - this note
   - the Chinese companion note
   - README links

The most important current conclusion is: the branch already demonstrates a
real SIMD-path win in the sort/preparation stage, but almost no win yet in the
actual CPAM tree build stage. That makes the next optimization target much
clearer.

## Practical Note

Generating a 1e9 dataset needs a large external data path. The repository's `.env.example` explicitly warns that the data directory should have `100GB+` free space.
