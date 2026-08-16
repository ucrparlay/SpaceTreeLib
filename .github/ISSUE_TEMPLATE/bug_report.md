---
name: Bug report
about: Something the library gets wrong
---

**What happened, and what you expected instead**

**A program that shows it**

The smallest one that still fails. Please say which tree, which split rule,
the point type and dimension, and roughly how many points.

**How it was built**

```
cmake -S . -B build ...      # the flags you used
g++ --version                # or clang
```

**Anything from the checks**

Does `ctest --test-dir build` pass? Does the failure survive
`-DPSI_ASSERTS=ON`, which keeps the assertions in an optimised build?
