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
`example`. It must never reach `include/parlaylib`, `include/libmorton` or
`include/highway` — those are submodules.

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

There is no automated correctness gate. To check a change:

- all nine targets build (eleven with `-DCGAL=ON`, which adds the
  `kd_ccp` / `p_ccp` oracles that `script/checkCorrect.sh` drives)
- the three `*_tree_example` binaries run and agree
- build with `-DDEBUG=ON` so the asserts and `validate()` run, then compare
  against a worktree at the previous commit on the same input

`orth_tree` requires a `spatial_median` partition rule. Pairing it with
`object_median` is now a `static_assert`, not a segfault.

`script/checkCorrect.sh` sets `tester` twice, so the second assignment wins
and `p_ccp` never runs: `p_tree` is not covered by it. Its `tag=0` also means
no batch update is exercised.

While iterating, check one size (`n=1000000`) — the larger sizes take over an
hour each in a `-DDEBUG=ON` build. Run the whole matrix once the work is
finished. Either way run `p_ccp` yourself, since the script will not.
