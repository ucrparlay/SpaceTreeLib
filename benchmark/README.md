# Recorded results

These 47 files are the summarized measurements behind the two papers' figures.
They are outputs, not inputs: the drivers under `tests/` write raw logs, the
`script_ae/merge_*.py` scripts read those logs and produce tables like these,
and the `plot_*.R` scripts turn the tables into figures.

| Directory | Input distribution |
|---|---|
| `uniform/`, `uniform_bigint/` | points uniform over the cube; `_bigint` when `-axis_max` exceeds 10^6 |
| `ss_varden/`, `ss_varden_bigint/` | skewed ("varden") distribution |
| `real_world/` | the real-world datasets |

Each row is one input file; the remaining columns are timings in seconds,
consumed positionally by the merge scripts. There is no header, and the column
meaning depends on which driver and which query mode produced it — read the
matching `script_ae/merge_*.py` to interpret a particular file.

## What is missing

**None of these files records the machine, the core count, the compiler, the
date, or the commit they were measured at**, and the repository has no tags,
so there is nothing to tie a number here to a version of the code. That
matters because the library is under active change: several recent commits
touch measured paths and say so, but there is no baseline to compare against.

If you know the provenance of these files, please record it here. For anything
measured from now on, the intent is:

- tag the commit a measurement campaign was run at, and name the tag here;
- record machine, core count, compiler and date alongside the numbers.

Timing runs must also be built with `-DPSI_NATIVE_ARCH=ON` (the default) on the
machine doing the measuring — see
[docs/ARTIFACT_EVALUATION.md](../docs/ARTIFACT_EVALUATION.md).
