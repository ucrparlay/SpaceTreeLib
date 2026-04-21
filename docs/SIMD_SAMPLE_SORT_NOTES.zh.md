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
2. 输出条目统一成紧凑的 `(curve_code, address)` 形式。
3. 如果 SIMD 路径没开，仍然可以退回原来的 CPAM sort 流程。

重点新增逻辑：

- `make_cpam_output_entry(...)`
  统一处理多种输出类型，包括原始 entry、augmentation key、传统 pair，以及新的整型 KVPair。
- `cpam_sample_sort_simd_pair(...)`
  现在不再先全局预计算整段 `code`。它会先准备 `tmp` 缓冲区，然后通过一个
  `LazyCurveCodePrepare` callback，在 SIMD sample sort 的顶层流程里按需：
  1. 绑定 payload 指针
  2. 计算 Morton code
  3. 回写排序 key
- `cpam_sample_sort_materialized_pair(...)`
  提供一个非 SIMD 的“先物化 code 再排序”的后备路径。
- [include/SIMD-Sample-Sort/src/two_pass/two_pass_simd.hpp](/home/xwang605/SpaceTreeLib/include/SIMD-Sample-Sort/src/two_pass/two_pass_simd.hpp:112)
  新增了 `prepare callback` 机制。这个 callback 只在 `Depth == 0` 生效：
  sample 阶段和 tile local sort 之前会触发一次，递归下去之后不再重复算 code。

从性能视角看，这次改动的关键已经不再是“先整段 precompute code，再 sort”，而是：

- legacy 路径：`fill code` 本来就是穿插在 sample sort 过程中做的
- 当前 SIMD 路径：也改成尽量贴近这种模式，在顶层 sample/tile 流程里懒绑定
  payload 并计算 code，而不是额外做一遍完整的大数组 `fill` pass

更准确地说，两条路线现在可以这样理解：

- `SIMD ON`
  是 `tmp` 先做 first-touch，然后在顶层 sample/tile 处理中懒执行
  `bind payload + compute code`，再走 SIMD sample sort，最后 build tree。
- `SIMD OFF`
  走的是原来的 CPAM sample sort 路线，`fill code` 不是一个完全独立的阶段，而是穿插在 sample sort 过程中完成的。

所以如果硬拆步骤：

- `SIMD ON`：更接近 `touch -> lazy(bind+code)+sort -> build`
- `SIMD OFF`：更接近 `fill+sort(interleaved) -> build`

这里要特别说明一点：当前 SIMD 路径之所以这样改，是为了更贴近原版
`legacy-interleaved-fill-sort` 的访存方式。如果先全局把所有 payload 绑定一遍，再 later
计算 code，通常会多出一遍完整的大数组写流量，不够像原版，也往往不利于性能。

## 新增的 Highway 2D materialize fast path

在当前这版代码里，`PTree` 的 2D benchmark 默认已经切到 `uint32_t x + uint32_t y + uint64_t code`
这套 `32/32/64` 布局。基于这个布局，`cpam_sample_sort` 额外加了一条仅针对 2D Morton 的
Highway fast path，位置在
[include/psi/dependence/cpam/cpam_sample_sort.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/cpam_sample_sort.h:247)。

这条路径只在下面这些条件同时满足时才会启用：

- `Point::Coord == uint32_t`
- `Point::GetDim() == 2`
- `Point::AT::CurveCode == uint64_t`
- `sizeof(Point) == 16`
- `sizeof(key_entry_pointer) == 16`

也就是说，它是一个很明确的“benchmark 场景专用 fast path”，不是对所有 `Point` 都生效的泛化实现。

实现思路是：

- 当前点数组在内存里就是交错的两列 `uint64_t word`
  - 第 0 列：打包后的 `(x, y)`，低 32 位是 `x`，高 32 位是 `y`
  - 第 1 列：`code`
- Highway 路径会在 `prepare_range(...)` 里按 lane 批量做：
  1. 取一批 packed `(x, y)`
  2. 用向量化的 `part1by1` / magic-bits 方式展开成 Morton code
  3. 回写 point 的 `code`
  4. 一次性写好 `tmp` 里的 `(key=code, payload=point_address)`

这里故意没有去改 `two_pass_simd.hpp` 的主排序流程，而是继续复用它已有的 `prepare callback` 机制：

- `sample` 阶段仍然走标量 `prepare_entry(...)`
- `tile local sort` 和顶层 base-case 这种区间准备阶段，才走 Highway 的 `prepare_range(...)`

这样做的原因是：

- 改动面最小
- 不会破坏现有 SIMD sample sort 的递归结构
- 很适合先验证“tile 级批量 materialize”本身到底能不能带来收益

如果触发了这条路径，`CPAM build` 的 route 会显示成：

```text
simd-pretouch-hwy2d-lazybindcode-sort
```

如果条件不满足，则仍然退回现有的：

```text
simd-pretouch-lazybindcode-sort
```

目前这里只实现了代码路径，还没有把新的 benchmark 数字写进文档；后续需要在新布局上重新测试。

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
- `PTree` 的 benchmark 入口现在默认使用 `uint32_t` 坐标
  - 2D 点的物理布局因此变成 `uint32_t x + uint32_t y + uint64_t code`
  - 也就是更贴近 Morton 实际使用的 `32/32/64` 模式，而不是之前的 `long/long/uint64_t`
- 为 `PTree` 输入增加按坐标去重逻辑，避免重复坐标影响构建
- 校验逻辑不再默认要求 `aug.id` 存在
- 新增 build-only 模式，并且现在专门走一个纯 build benchmark helper

这个改动放在 [tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:1899) 和
[tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:2061)：

- benchmark 读点时，`unsigned` 坐标会走显式范围检查后再转换
- `PTree` 的 2D/3D benchmark 点类型改成了 `AugPoint<uint32_t, dim, AugCode>`

这样 legacy 路线、当前 SIMD 路线，以及后续要加的 Highway materialize 路线，都会基于同一套更紧凑、也更符合 Morton encode 语义的点布局做比较。

你现在加的 build-only 逻辑在 [tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:1245)：

```cpp
if (kTag == 0 && kQueryType == 0) {
  BenchmarkBuildOnly<Point, Tree>(wp, tree, kRounds);
  return;
}
```

也就是说，现在只要传 `-t 0 -q 0`，框架就只会建树，不再进入 update 或 query 路径。

这里最关键的一点是：build-only 不再复用通用的 `BuildTree(...)`。因为
`BuildTree(...)` 原本是给后续 insert/query 流程服务的，计时结束后还会额外
再建一次树，把状态留给后续步骤。对于纯 build benchmark，这会让日志里多出
一条“不属于统计口径本身”的 build。

现在 `BenchmarkBuildOnly(...)` 的语义是：

1. 先做 1 次 warmup `build + delete`
2. 再做 3 次正式的 `build + delete`
3. 只统计后面 3 次的平均值

这样就和“先丢掉第一轮，再平均后三轮”的预期一致了。

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

如果你想切换 Morton encode 的实现路径，当前也支持一个单独的 CMake 开关：

```bash
cmake -S . -B build-simd \
  -DUSE_SIMD_SAMPLE_SORT=ON \
  -DUSE_LIBMORTON_SIMD_ENCODE=ON
```

其中：

- `USE_LIBMORTON_SIMD_ENCODE=ON`
  使用 libmorton 的默认 encode 入口
- `USE_LIBMORTON_SIMD_ENCODE=OFF`
  退回之前显式指定的 `m2D_e_for / m3D_e_for_ET`

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
[CPAM build] route=simd-pretouch-lazybindcode-sort n=... touch=... fill=0.000000 sort=... build=... total=...
```

或者：

```text
[CPAM build] route=legacy-interleaved-fill-sort n=... fill=interleaved-with-sort sort+fill=... build=... total=...
```

含义是：

- `route`
  标明当前走的是哪条 sort/build 路线
- `touch`
  只在新的 SIMD 路线下单独打印，表示 `tmp` 的并行 first-touch / 预热时间
- `fill`
  现在更多表示“独立预处理 pass”的时间。对当前 `lazybindcode` 路线来说，这个值通常是 `0`
  或接近 `0`，因为 payload 绑定和 code 计算已经并入 sort 顶层流程了。
- `sort`
  对当前 SIMD lazy 路线来说，这里面除了纯排序内核，也包含顶层 sample/tile 阶段内的
  `bind payload + compute code`
- `sort+fill`
  旧路线里 fill 和 sort 无法完全拆开，所以合并打印
- `build`
  从排好序的结果继续建树的时间

## 当前 1e9 Morton 的实测结果

先放当前代码，也就是 `simd-pretouch-lazybindcode-sort` 这版的最新结果。

测试命令：

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

口径说明：

- 第 1 条 `[CPAM build]` 是 warmup
- 后 3 条 `[CPAM build]` 是正式统计轮次
- `p_test` 最后一行打印的标量，是 build-only harness 对正式轮次求出来的平均值

对后 3 轮求平均后，结果如下：

| 路线 | fill / sort 阶段 | build 阶段 | CPAM total | build-only 平均值 |
| --- | ---: | ---: | ---: | ---: |
| 原版 (`legacy-interleaved-fill-sort`) | `1.542110 + 1.568069 + 1.534770` 平均 `1.548316s` (`sort+fill`) | `0.680652 + 0.680876 + 0.679833` 平均 `0.680454s` | `2.228770s` | `2.562340s` |
| SIMD (`simd-pretouch-lazybindcode-sort`) | `0.024810s touch + 1.346630s sort` | `0.680346s` | `2.051785s` | `2.436810s` |

这里 SIMD 路线里的 `touch` 和 `sort` 分别是：

- `touch`
  `(0.025256 + 0.024774 + 0.024400) / 3 = 0.024810s`
- `sort`
  `(1.423008 + 1.310194 + 1.306687) / 3 = 1.346630s`
- `build`
  `(0.679280 + 0.680391 + 0.681366) / 3 = 0.680346s`
- `CPAM total`
  `0.024810 + 1.346630 + 0.680346 = 2.051785s`

差异可以总结成：

- `[CPAM build] total`
  从 `2.228770s` 降到 `2.051785s`，快了 `0.176985s`
  (`1.086x`，提升约 `7.94%`)
- `build-only` 最后一行标量
  从 `2.562340s` 降到 `2.436810s`，快了 `0.125530s`
  (`1.052x`，提升约 `4.90%`)

### 为什么 `[CPAM build]` 能快 0.1 秒，但最后标量只快 0.04 秒

这个不是计时器坏了，而是**两种打印口径不一样**。

`[CPAM build]` 这行是在 [map.h](/home/xwang605/SpaceTreeLib/include/psi/dependence/cpam/map.h:364)
里面打的，它覆盖的是：

1. `Build::sort_remove_duplicates(SS)` 里的 sort 路径
2. `Tree::multi_insert_sorted(...)` 的 CPAM build 路径

也就是说，它主要反映的是 **CPAM 内核这段**。

但 `p_test` 最后一行的标量来自 [tests/test_framework.h](/home/xwang605/SpaceTreeLib/tests/test_framework.h:366)：

```cpp
t.start();
pkd.Build(wp.cut(0, n));
total_build += t.next_time();
```

这里计的是**整个 `pkd.Build(...)`**。对 `PTree` 来说，`Build(...)` 外面还有一层：
[p_build_tree.hpp](/home/xwang605/SpaceTreeLib/include/psi/p_tree_impl/p_build_tree.hpp:23)

```cpp
auto aux = Points::uninitialized(parlay::size(In));
parlay::copy(In, parlay::make_slice(aux));
Slice A = parlay::make_slice(aux);
Build_(A);
```

也就是说，最后那个 build-only 标量还额外包含：

- 一次 `aux` 大数组分配
- 一次把输入整体拷贝到 `aux` 的 full copy

这部分开销在 SIMD 和 legacy 两条 CPAM 路线上基本都要付，而且不在 `[CPAM build]`
打印里，所以会把 CPAM 内核里的收益“稀释”掉。

从你这组数据也能直接看出来：

- SIMD: `2.436810 - 2.051785 = 0.385025s`
- Legacy: `2.562340 - 2.228770 = 0.333570s`

这大约 `0.33s ~ 0.41s` 的区间，就是 `PTree::Build(...)` 外层包装开销和一些运行时噪声。
所以：

- 如果你想看“SIMD sort / CPAM build 核心路径到底快了多少”，看 `[CPAM build]`
- 如果你想看“用户调用一次 `tree.Build(...)` 端到端到底快了多少”，看最后那个标量

下面这组表格是更早一版 `simd-precompute-fill-sort` 路线的历史结果，保留它主要是为了记录
“先分离 fill/sort” 时的旧口径。

| 路线 | fill / sort 阶段 | build 阶段 | CPAM total |
| --- | ---: | ---: | ---: |
| 原版 (`legacy-interleaved-fill-sort`) | `1.811655s` (`sort+fill`) | `0.692736s` | `2.504391s` |
| SIMD (`simd-precompute-fill-sort`) | `0.356934s fill + 1.145278s sort = 1.502212s` | `0.685634s` | `2.187846s` |

差异可以直接总结成：

- `fill+sort` 快了 `0.309443s`，从 `1.811655s` 降到 `1.502212s`
  (`1.206x`，提升 `17.08%`)
- `build` 只快了 `0.007102s`，从 `0.692736s` 到 `0.685634s`
  (约 `1.03%`，基本可以视为持平)
- CPAM `total` 快了 `0.316546s`，从 `2.504391s` 到 `2.187846s`
  (`1.145x`，提升 `12.64%`)

所以现在这条 SIMD 路线的优势，几乎都集中在建树前的准备阶段，也就是
`fill + materialize + sort`；真正的 tree build kernel 本身，目前几乎没变快。

如果看 `p_test` 最后一行的 build-only 平均值：

- 原版：`2.85636`
- SIMD：`2.63982`

也就是端到端的 build-only 平均值快了 `0.21654s`
(`1.082x`，提升 `7.58%`)。这个提升比 `[CPAM build]` 里的 `total` 更小，是因为
它还混入了 trace 之外的 framework 开销。

## 当前分支改动总览

到目前为止，这个 branch 的改动可以归成 5 组：

1. CPAM sort / data-path 改动
   - `PTree` 改成物化 `(CurveCode, address)` 形式
   - 新增 SIMD sample sort 适配层
   - 新增非 SIMD 的 materialized fallback
   - 新增 route / timing trace
2. CPAM build 兼容性改动
   - 让构建阶段同时支持传统指针 payload 和整数化地址 payload
3. 测试框架改动
   - `AugCode`
   - 坐标去重
   - 不再依赖 `aug.id` 的校验逻辑
   - 独立的 build-only benchmark 流程
4. 构建与集成改动
   - `USE_SIMD_SAMPLE_SORT`
   - `PRINT_CPAM_BUILD_TIMING`
   - `include/SIMD-Sample-Sort` / Highway 集成
   - `.gitignore` 忽略 `build-simd/`
5. 文档改动
   - 这份说明
   - 英文说明
   - README 链接

目前最重要的结论是：这个分支已经证明 SIMD 路线在“排序/准备阶段”有真实收益，
但在“实际 CPAM build 阶段”还几乎没有收益。所以下一步优化目标已经很明确了。

## 实际使用提醒

仓库里的 `.env.example` 明确写了数据目录需要 `100GB+` 可用空间。所以 1e9 数据建议直接生成到外部数据盘，不适合放在当前仓库目录或 `/tmp` 这种空间有限的位置。
