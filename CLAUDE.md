# PSI

Header-only C++20 parallel spatial indexes: kd_tree, orth_tree, p_tree.
Backs two papers (SIGMOD'25, PPoPP'26), so measured numbers matter.

## Coding style

Kernel formatting, C++ standard library naming. When the two disagree with
the existing code, match the file you are editing and say so.

Formatting is whatever `.clang-format` produces: tabs, 8 wide, 80 columns,
Linux braces. Run `make format` and never argue with it. Every file must be a
clang-format fixpoint before commit.

`make format` covers `include/psi` (including `dependence/cpam`), `tests` and
`example`. It must never reach `include/parlaylib` or `include/libmorton` —
those are submodules.

Comments are `/* */`. The `//` still in the tree is legacy.

Write the simple version. No abstraction until there are three callers.
No helper that is used once. No option nobody asked for.

## Naming

As the standard library names things.

- classes, functions, variables, concepts: `snake_case`
- member type aliases: `snake_case` ending in `_type` (`box_type`, `slice_type`)
- constants: `snake_case`, no prefix (`bucket_num`, not `kBucketNum`)
- template parameters: `CamelCase` (`Point`, `SplitRule`, `SkHeight`)
- private data members: trailing underscore (`root_`)
- macros: `UPPER_CASE`

A few names stay `CamelCase` because a third party owns them: `Tree` is
CPAM's and CGAL's nested type, and `kd_ccp.cpp` typedefs CGAL names onto
themselves (`Point_d`, `Fuzzy_iso_box`). Renaming those breaks the build in
ways the compiler reports far from the cause.

## Configuration

A tree takes one template parameter, a traits type:

    using traits = psi::tree_traits<point, split_rule>;
    psi::kd_tree<traits> tree;

`psi/tree_traits.h` holds both halves. `point_traits<Point>` is what follows
from the point type alone -- the aliases, the tuning constants, and the box
arithmetic; the default augmentations take it, so it cannot name the traits
that name them. `tree_traits<Point, SplitRule, LeafAug, InteriorAug,
SkHeight, ImbaRatio>` adds the rest and inherits it. `orth_tree_traits` is
the same thing with a skeleton height that fits the dimension.

`base_tree` inherits its traits, so a new alias or constant needs adding in
one place. A dependent base is invisible to unqualified lookup, so anything
a tree uses is re-exported by the block at the top of `base_tree.h`; a name
missing from that block is a compile error, not a silent fallback.

## Comments

Few and short. A comment costs maintenance; most code does not earn one.

Comment why, not what. If the code needs a comment to say what it does,
rename something instead.

    /* Points on the split plane go right. */

Do not narrate the diff. Do not write "this used to be X" or restate the
commit message. Git remembers.

No decoration: no `NOTE:`, `WARN:`, `PERF:` tags, no banner lines, no
comment blocks that repeat the function signature.

Delete commented-out code. Git remembers that too.

## Words

Plain English. Say the thing.

Do not invent terms and do not bend a word away from its meaning. If a
function marks a node for rebuilding, say marks, not "tombs". If it
partitions points, say partitions.

Existing bent names (`Seieve`, `Sparcy`, `tomb`, `puffy`) are legacy. Do not
add more. Renaming them is fine as its own change.

No jargon and no filler: not "leverage", "robust", "seamless", "simply".

## Do not touch

The parallel algorithms are frozen. Changing their control flow, work, or
depth is out of scope unless asked:

- the eight `*_recursive` functions in `kd_tree_impl/` and `orth_tree_impl/`
- the skeleton tag values `bucket_num+1/+2/+3` in `inner_tree.hpp`
- the force-parallel flag and the `serial_build_cutoff` fallback
- the granularity argument on every `parallel_for`
- `-march=native`, `-O3`, `-mcx16`

`include/baselines/` is out of scope. It exists to compare against; do not
refactor it. Only its two PSI-facing files (`cpam_raw/cpamtree.hpp`,
`cpam_raw/cpam/augmented_ops.h`) follow our names, and only so they compile.

`include/psi/dependence/cpam/` is a permanent PSI fork of CPAM, not vendored
code — it diverges from upstream by ~2,800 lines, and upstream is no longer
maintained, so there is no re-sync to preserve changes for. It is nonetheless
deferred: leave it alone unless asked. See `THIRD_PARTY.md`.

Two consequences of that deferral, so they are not rediscovered:

- The 7 `-Wmaybe-uninitialized` warnings in a release build all come from
  `dependence/cpam/augmented_node.h` (a partly indeterminate `std::pair`).
  They are benign, but they will show up in any sanitizer run touching
  `p_tree`. "`include/psi` is warning-clean" means *apart from these*.
- `p_tree` has no `validate()`; its node type is CPAM's, which the validation
  branches do not cover. Cover `p_tree` with black-box tests instead.

Flag anything that could move a measured number, even an improvement.

## Build and check

    cmake -S . -B build && make -C build -j

Warnings are on for the three examples only; they cover the whole library.
`include/psi` must stay warning-clean apart from the vendored
`space_filling_curve/hilbert.c`.

The tree is formatted with **clang-format 19.1.0**, which is what CI checks
against. Any version from 20 on rewrites `space_filling_curve/hilbert.c`
wholesale, and it is not a disagreement about style: clang-format 20 made `C`
a language kind separate from `Cpp` and infers it from the extension, and that
file is C++ -- `namespace psi`, `inline` functions, `#include`d by
`splitter.h`. Parsed as C, `namespace psi` is an identifier before a brace
block, so the body gets indented. The same bytes named `.hpp` are already a
fixpoint under 22. Until the file is renamed, format with 19.1.0.

There is a correctness gate now: `ctest` in the build directory runs the unit
suite under `tests/unit`, which compares all three trees against brute force
(queries, batch updates, flatten round-trip, degenerate inputs, double
coordinates). It takes well under a second and needs no data files. Unit tests
include only `psi/…` headers -- never `tests/test_framework.h`, which costs
~70 s and 1.7 GB per translation unit and pulls in the baselines.

`-DPSI_ASSERTS=ON` keeps `assert()` in an optimised build, which is the
configuration to check a change under; `-DDEBUG=ON` also works but is `-O0`.

To check a change:

- `ctest` passes
- all nine targets build (eleven with `-DCGAL=ON`, which adds the
  `kd_ccp` / `p_ccp` oracles that `script/checkCorrect.sh` drives)
- the three `*_tree_example` binaries run and agree -- they check their own
  answers against brute force now, so a non-zero exit means a wrong answer
- build with `-DDEBUG=ON` so the asserts and `validate()` run, then compare
  against a worktree at the previous commit on the same input

`orth_tree` requires a `spatial_median` partition rule. Pairing it with
`object_median` is now a `static_assert`, not a segfault.

`script/checkCorrect.sh` compares against CGAL and exits non-zero on the first
mismatch: `./checkCorrect.sh [fast|full] [data-root] [build-dir]`. It sweeps
all three trees, the update bits and the query bits. Generate your own data
with `data_generator` and point it at the root, or set `PSI_DATA`; `PSI_NODES`
overrides the sizes. Note `-t` bit 1 (delete) needs bit 0 (insert): alone it
would ask `batch_delete` to remove points that were never inserted, which its
contract forbids.

While iterating, check one size (`n=1000000`) — the larger sizes take over an
hour each. Run the whole matrix once the work is finished.

`data_generator` writes `uniform_bigint/` and `ss_varden_bigint/` whenever the
axis maximum exceeds 10^6, but the script looks for `uniform/` and
`ss_varden/`; symlink them or the sweep finds no files and exits non-zero.
