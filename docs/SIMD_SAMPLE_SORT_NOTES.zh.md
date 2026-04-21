# SIMD Sample Sort 分支说明

这份文档总结了 `sean/sort` 分支上的主要改动、开启 SIMD sample sort 的编译方式，以及当前只跑建树的测试命令。

## 这次改动的核心

这次分支最核心的变化，是 `PTree` 在构建阶段不再沿用旧的 `(AugType, Point*)` 排序形式，而是先把待排序项物化成更紧凑的键值结构：

- key：`CurveCode`，类型是 `uint64_t`
- value：点的地址，存成 `uint64_t`

对应的类型切换在 [include/psi/p_tree.h](/home/xwang605/SpaceTreeLib/include/psi/p_tree.h:66)：

```cpp
using key_entry_pointer =
    parlay::KVPair<typename Point::AT::CurveCode, uint64_t>;
```

这样做之后，CPAM 排序面对的就是“纯整数 key + 纯整数 payload”，和 SIMD sample sort 的接口更匹配，也更适合走向量化排序实现。

## sort 路径是怎么改的

[include/psi/dependence/cpam/cpam_sample_sort.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/cpam_sample_sort.h:33) 现在新增了几层适配：

1. 通过 `PSI_USE_SIMD_SAMPLE_SORT` 判断是否启用 SIMD 排序。
2. 排序前先预计算 filling-curve code，并生成 `(code, address)` 数组。
3. 如果 SIMD 路径没开，仍然可以退回原来的 CPAM sort 流程。

重点新增逻辑：

- `make_cpam_output_entry(...)`
  统一处理多种输出类型，包括原始 entry、augmentation key、传统 pair，以及新的整型 KVPair。
- `cpam_sample_sort_simd_pair(...)`
  先物化 `(curve_code, address)`，再调用 SIMD sample sort。
- `cpam_sample_sort_materialized_pair(...)`
  提供一个非 SIMD 的“先物化 code 再排序”的后备路径。

从性能视角看，这次改动本质上是把原来“排序时不断处理复杂对象”的过程，改成“先把可比较的整数 key 准备好，再对紧凑记录排序”。

更准确地说，两条路线现在可以这样理解：

- `SIMD ON`
  是先单独做 `fill code + materialize KV`，然后单独做 SIMD sort，最后再 build tree。
- `SIMD OFF`
  走的是原来的 CPAM sample sort 路线，`fill code` 不是一个完全独立的阶段，而是穿插在 sample sort 过程中完成的。

所以如果硬拆步骤：

- `SIMD ON`：`fill -> sort -> build`
- `SIMD OFF`：更接近 `fill+sort(interleaved) -> build`

## 为什么 CPAM build 也要一起改

因为现在 payload 不再直接是 `Point*`，而是整数化后的地址，所以 CPAM build 里凡是“从排序结果里恢复原始 entry”的地方都要一起适配。

相关文件有：

- [include/psi/dependence/cpam/basic_node_helpers.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/basic_node_helpers.h:11)
- [include/psi/dependence/cpam/sequence_ops.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/sequence_ops.h:721)
- [include/psi/dependence/cpam/build.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/build.h:42)

这几处的改法都比较一致：通过统一的 helper，把“传统指针 payload”和“整型地址 payload”两种形式都还原成原始 entry，保证后续建树逻辑不用大改。

## 测试框架的配套改动

[tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:73) 里也做了几项重要配套：

- 新增 `AugCode`，只保存 curve code，不再依赖 `id`
- 为 `PTree` 输入增加按坐标去重逻辑，避免重复坐标影响构建
- 校验逻辑不再默认要求 `aug.id` 存在
- 新增 build-only 模式

你现在加的 build-only 逻辑在 [tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:1245)：

```cpp
if (kTag == 0 && kQueryType == 0) {
  BuildTree<Point, Tree, kTestTime>(wp, kRounds, tree);
  tree.DeleteTree();
  return;
}
```

也就是说，现在只要传 `-t 0 -q 0`，框架就只会建树，不再进入 update 或 query 路径。

## 开启 SIMD sort 的编译命令

在子模块准备好的前提下，推荐直接这样编译：

```bash
cmake -S . -B build-simd -DUSE_SIMD_SAMPLE_SORT=ON
cmake --build build-simd -j
```

如果你还想打开 CPAM build 的阶段计时打印，可以一起开：

```bash
cmake -S . -B build-simd \
  -DUSE_SIMD_SAMPLE_SORT=ON \
  -DPRINT_CPAM_BUILD_TIMING=ON
cmake --build build-simd -j --target p_test data_generator
```

如果你只关心这次要用到的两个可执行文件，也可以只编译它们：

```bash
cmake --build build-simd -j --target p_test data_generator
```

如果你想在编译时启用 jemalloc，对应的 CMake 选项是：

```bash
cmake -S . -B build-simd \
  -DUSE_SIMD_SAMPLE_SORT=ON \
  -DJEMA=ON
```

这里不是 `jlloc`，仓库里这个开关名字叫 `JEMA`，对应的是 jemalloc。

## 生成 1e9 的 uniform 数据

生成器文档在 [docs/MANUAL.md](/home/xwang605/SpaceTreeLib/docs/MANUAL.md:232)，实现代码在 [tests/data_generate.cpp](/home/xwang605/SpaceTreeLib/tests/data_generate.cpp:1)。

如果你只想生成一份 1e9 的 uniform 数据，不做 data washing / sort，可以直接：

```bash
./build-simd/data_generator \
  -p /path/to/data \
  -n 1000000000 \
  -d 2 \
  -file_num 1 \
  -varden 0
```

生成结果会落在：

```text
/path/to/data/uniform_bigint/1000000000_2/1.in
```

如果你跑的是 3D，并且希望和现有 AE 脚本一致，用较小坐标范围：

```bash
./build-simd/data_generator \
  -p /path/to/data \
  -n 1000000000 \
  -d 3 \
  -file_num 1 \
  -varden 0 \
  -axis_max 1000000
```

## 当前 build-only 的运行命令

你现在这个 framework 改动下，SPaC-tree 的 build-only 命令可以写成：

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

参数说明：

- `-T 2`：选择 `PTree`
- `-l 2`：Morton curve
- `-l 1`：Hilbert curve
- `-t 0 -q 0`：只建树
- `-i 0`：不读取插入数据

## 计时打印格式

打开 `-DPRINT_CPAM_BUILD_TIMING=ON` 之后，从空树建 `PTree` 时会额外打印一行：

```text
[CPAM build] route=simd-precompute-fill-sort n=... fill=... sort=... build=... total=...
```

或者：

```text
[CPAM build] route=legacy-interleaved-fill-sort n=... fill=interleaved-with-sort sort+fill=... build=... total=...
```

含义是：

- `route`
  标明当前走的是哪条 sort/build 路线
- `fill`
  只在 SIMD / precompute 路线下单独存在
- `sort`
  纯排序时间
- `sort+fill`
  旧路线里 fill 和 sort 无法完全拆开，所以合并打印
- `build`
  从排好序的结果继续建树的时间

## 实际使用提醒

仓库里的 `.env.example` 明确写了数据目录需要 `100GB+` 可用空间。所以 1e9 数据建议直接生成到外部数据盘，不适合放在当前仓库目录或 `/tmp` 这种空间有限的位置。
