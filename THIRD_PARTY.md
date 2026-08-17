# Third-party code in PSI

The root `LICENSE` (MIT) covers PSI's own source. It does not cover the
components below, which are redistributed under their own terms. This file
exists so a license scan, or a user shipping a product built on PSI, can see
what is here and what it requires.

## Required notice

`include/psi/dependence/space_filling_curve/hilbert.hpp` is copyrighted by Rice
University. Its terms require the following to be displayed anywhere the
end-user would see PSI's own copyright:

> Hilbert Curve implementation copyright 1998, Rice University

The same terms ask that modifications carry the name of the person who made
them, the date, and the reason. Recorded here rather than in the file, so its
notice stays verbatim:

- made `inline`, so two translation units can both include it and still link;
- wrapped in `namespace psi::hilbert`;
- renamed from `hilbert.c` to `hilbert.hpp` on 2026-08-17. It was never C --
  the namespace and the `inline` see to that -- and clang-format 20 began
  inferring `Language: C` from the extension and rewriting the whole file.
  Its own header still names it `hilbert.c`, which is the upstream file.

## Git submodules

Not distributed in this repository — fetched by `git submodule update --init`.

| Path | Upstream | Pinned | License |
|---|---|---|---|
| `include/parlaylib` | https://github.com/cmuparlay/parlaylib | `36459f4` (master) | MIT |
| `include/libmorton` | https://github.com/Forceflow/libmorton | `7923faa` (v0.2.12) | MIT |

## Vendored copies

These ship in the repository.

### `include/psi/dependence/cpam/` — a permanent fork of CPAM

Derived from CPAM (https://github.com/cmuparlay/CPAM, MIT), which builds on
PAM (https://github.com/cmuparlay/PAM, MIT). **This is a fork, not a mirror.**

Measured against the pristine copy under `include/baselines/cpam_raw`
(whitespace-insensitive, at commit `da195e9` so PSI's later renaming is
excluded): it diverges by roughly **2,842 lines** across the shared files,
with `augmented_ops.h` rewritten in 481 of 562 lines, `map_ops.h` in 882 and
`map.h` in 460. It also adds `cpam_sample_sort.h`, which applies the
space-filling curve. Those changes are PSI's own work and carry PSI's
copyright.

Upstream CPAM/PAM is no longer maintained, so this fork is permanent. Do not
plan to re-sync it, and do not read the divergence as drift to be reconciled.

The exact upstream commit these files were taken from is not recorded anywhere
in the tree and could not be recovered; if you know it, add it here.

### `include/baselines/cpam_raw/` — CPAM and PAM, kept pristine

A second, unmodified copy of the same upstream projects (MIT), retained as the
benchmark baseline that PSI is measured against. Only two files were touched,
and only so they keep compiling against PSI's renamed types: `cpamtree.hpp`
and `cpam/augmented_ops.h`.

Note for anyone comparing the two copies: PSI's fork reads a node's augmented
bounding box by `const &` in the knn descent (`augmented_node.h`, used from
`augmented_ops.h`) where this pristine copy returns it by value. That
difference sits on the path the SPaC-vs-CPAM figures measure. It predates the
2026 refactor and is left as-is; it is recorded here so it is not mistaken for
an accident.

### `include/baselines/zdtree/`, `include/baselines/zdtree_3d/`

Implementations of a competing design, used only for benchmarking. They derive
from PSI's own `base_tree`, and no upstream source or third-party copyright
notice was found in them. If any of this code came from elsewhere, it should be
attributed here.

## Scope note

Everything under `include/baselines/` exists to compare against and is not
maintained as part of the library.
