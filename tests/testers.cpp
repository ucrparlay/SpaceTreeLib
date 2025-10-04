#include <tuple>

#include "test_framework.h"

template <typename Tree>
struct InnerNode {
  using Box = typename Tree::Box;
  InnerNode()
      : box(Tree::GetEmptyBox()), leaf_offset(0), size(0), is_leaf(false) {}
  InnerNode(Box const& _box, size_t _leaf_offset, size_t _size, bool _is_leaf)
      : box(_box), leaf_offset(_leaf_offset), size(_size), is_leaf(_is_leaf) {}

  Box box;
  size_t leaf_offset;
  size_t size;
  bool is_leaf;
};

template <typename Tree>
void LinearDecomposition(Node* T, auto& leaf_seq, auto& inner_tree_seq,
                         size_t& leaf_offset, size_t inter_idx) {
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  if (T->is_leaf) {
    Leaf* TL = static_cast<Leaf*>(T);
    assert(!TL->is_dummy);
    inner_tree_seq[inter_idx] = {TL->GetBox(), leaf_offset, T->size,
                                 T->is_leaf};

    for (int i = 0; i < T->size; i++) {
      leaf_seq[leaf_offset++] = static_cast<Leaf*>(T)->pts[i];
    }
    return;
  }
  Interior* TI = static_cast<Interior*>(T);
  inner_tree_seq[inter_idx] = {TI->GetBox(), leaf_offset, T->size, T->is_leaf};

  assert(inter_idx * 2 + 1 < inner_tree_seq.size());
  LinearDecomposition<Tree>(TI->left, leaf_seq, inner_tree_seq, leaf_offset,
                            inter_idx * 2);
  LinearDecomposition<Tree>(TI->right, leaf_seq, inner_tree_seq, leaf_offset,
                            inter_idx * 2 + 1);
  return;
}

template <typename Tree, typename Range, typename Box>
void RangeQueryLinear(size_t idx, auto const& leaf_seq,
                      auto const& inner_tree_seq, Range Out, size_t& s,
                      Box const& query_box, RangeQueryLogger& logger) {
  if (inner_tree_seq[idx].is_leaf) {
    auto const& node = inner_tree_seq[idx];
    for (int i = 0; i < node.size; i++) {
      if (Tree::WithinBox(leaf_seq[node.leaf_offset + i], query_box)) {
        Out[s++] = leaf_seq[node.leaf_offset + i];
      }
    }
    logger.vis_leaf_num++;
    return;
  }
  logger.vis_interior_num++;

  auto recurse = [&](size_t next_idx) -> void {
    auto next_node = inner_tree_seq[next_idx];
    auto const& box = next_node.box;
    if (!Tree::BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      return;
    } else if (Tree::WithinBox(box, query_box)) {
      logger.full_box_num++;
      // FlattenRec<Leaf, Interior>(Ts, Out.cut(s, s + Ts->size));
      for (int i = 0; i < next_node.size; i++) {
        Out[s++] = leaf_seq[next_node.leaf_offset + i];
      }
      return;
    } else {
      RangeQueryLinear<Tree>(next_idx, leaf_seq, inner_tree_seq, Out, s,
                             query_box, logger);
      return;
    }
  };

  recurse(idx * 2);
  recurse(idx * 2 + 1);
  return;
}

template <typename Point, typename Tree, typename Range>
void KNNLinear(size_t idx, auto& leaf_seq, auto& inner_tree_seq, Point const& q,
               kBoundedQueue<Point, Range>& bq, KNNLogger& logger) {
  using Num = typename Tree::Num;
  auto const& node = inner_tree_seq[idx];
  if (bq.size() &&
      Num::Gt(Tree::P2BMinDistanceSquare(q, node.box), bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }

  if (node.is_leaf) {
    logger.vis_leaf_num++;

    size_t i = 0;
    while (!bq.full() && i < node.size) {
      bq.insert(std::make_pair(
          std::ref(leaf_seq[node.leaf_offset + i]),
          Tree::P2PDistanceSquare(q, leaf_seq[node.leaf_offset + i])));
      i++;
    }
    while (i < node.size) {
      auto r = Tree::InterruptibleDistance(q, leaf_seq[node.leaf_offset + i],
                                           bq.top_value());
      if (Num::Lt(r, bq.top_value())) {  // PERF: the queue is full, no need to
                                         // insert points with equal distances
        bq.insert(std::make_pair(std::ref(leaf_seq[node.leaf_offset + i]), r));
      }
      // else if (TL->is_dummy) {
      //   break;
      // }
      i++;
    }
    return;
  }

  logger.vis_interior_num++;
  // Interior* TI = static_cast<Interior*>(T);
  Coord dist_left = Tree::P2BMinDistanceSquare(q, inner_tree_seq[idx * 2].box);
  Coord dist_right =
      Tree::P2BMinDistanceSquare(q, inner_tree_seq[idx * 2 + 1].box);
  bool go_left = Num::Leq(dist_left, dist_right);

  KNNLinear<Point, Tree, Range>(go_left ? idx * 2 : idx * 2 + 1, leaf_seq,
                                inner_tree_seq, q, bq, logger);

  logger.check_box_num++;
  if (Num::Gt(go_left ? dist_right : dist_left, bq.top_value()) && bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNLinear<Point, Tree, Range>(go_left ? idx * 2 + 1 : idx * 2, leaf_seq,
                                inner_tree_seq, q, bq, logger);
  return;
}

template <typename Point, typename Tree, typename Range>
auto KNNLinearWrapper(auto& leaf_seq, auto& inner_tree_seq, Point const& q,
                      kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  KNNLinear<Point, Tree, Range>(1, leaf_seq, inner_tree_seq, q, bq, logger);
  return logger;
}

template <typename Point, typename Tree, bool printHeight = 0,
          bool printVisNode = 1>
void QueryKNNLinear(auto& leaf_seq, auto& inner_tree_seq,
                    uint_fast8_t const& Dim, parlay::sequence<Point> const& WP,
                    int const& rounds, Typename* kdknn, int const K) {
  using Points = typename Tree::Points;
  using Coord = typename Point::Coord;
  using DisType = typename Point::DisType;
  using nn_pair = std::pair<std::reference_wrapper<Point>, DisType>;
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  // using nn_pair = std::pair<Point, DisType>;
  size_t n = WP.size();
  // int LEAVE_WRAP = 32;
  double loopLate = rounds > 1 ? 0.01 : -0.1;

  Points wp = WP;
  // Points wp = parlay::random_shuffle(WP);

  parlay::sequence<nn_pair> Out(
      K * n, nn_pair(std::ref(wp[0]), static_cast<DisType>(0)));
  // parlay::sequence<nn_pair> Out(K * n);
  parlay::sequence<kBoundedQueue<Point, nn_pair>> bq =
      parlay::sequence<kBoundedQueue<Point, nn_pair>>::uninitialized(n);
  parlay::parallel_for(
      0, n, [&](size_t i) { bq[i].resize(Out.cut(i * K, i * K + K)); });
  parlay::sequence<size_t> vis_leaf(n), vis_inter(n), gen_box(n), check_box(n),
      skip_box(n);

  double aveQuery = time_loop(
      rounds, loopLate,
      [&]() { parlay::parallel_for(0, n, [&](size_t i) { bq[i].reset(); }); },
      [&]() {
        parlay::parallel_for(0, n, [&](size_t i) {
          // for (size_t i = 0; i < n; i++) {
          auto [vis_leaf_num, vis_inter_num, gen_box_num, check_box_num,
                skip_box_num] =
              KNNLinearWrapper<Point, Tree>(leaf_seq, inner_tree_seq, wp[i],
                                            bq[i]);
          kdknn[i] = bq[i].top().second;
          vis_leaf[i] = vis_leaf_num;
          vis_inter[i] = vis_inter_num;
          gen_box[i] = gen_box_num;
          check_box[i] = check_box_num;
          skip_box[i] = skip_box_num;
          // }
        });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;
  if (printVisNode) {
    std::cout << static_cast<double>(parlay::reduce(vis_leaf.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(vis_inter.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(gen_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(check_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
    std::cout << static_cast<double>(parlay::reduce(skip_box.cut(0, n))) /
                     static_cast<double>(n)
              << " " << std::flush;
  }

  return;
}

template <typename Tree, typename Range, typename Box>
auto RangeQueryLinearWrapper(auto const& inner_tree_seq, auto const& leaf_seq,
                             Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  RangeQueryLinear<Tree>(1, leaf_seq, inner_tree_seq, parlay::make_slice(Out),
                         s, query_box, logger);
  return std::make_pair(s, logger);
}

template <typename Point, typename Tree>
void RangeQueryFixLinear(auto& leaf_seq, auto& inner_tree_seq, Typename* kdknn,
                         int const& rounds, parlay::sequence<Point>& Out,
                         int rec_type, int rec_num, int DIM,
                         auto const& query_box_seq, auto max_size) {
  parlay::sequence<size_t> vis_leaf(rec_num), vis_inter(rec_num),
      gen_box(rec_num), full_box(rec_num), skip_box(rec_num);
  auto [offset, tot_size] = parlay::scan(
      parlay::delayed_tabulate(
          rec_num, [&](size_t i) -> size_t { return query_box_seq[i].second; }),
      parlay::addm<size_t>());
  offset.push_back(tot_size);
  Out.resize(tot_size);

  double aveQuery = time_loop(
      rounds, 0.01, [&]() {},
      [&]() {
        parlay::parallel_for(0, rec_num, [&](size_t i) {
          auto [size, logger] = RangeQueryLinearWrapper<Tree>(
              inner_tree_seq, leaf_seq, query_box_seq[i].first,
              Out.cut(offset[i], offset[i + 1]));

          kdknn[i] = size;
          vis_leaf[i] = logger.vis_leaf_num;
          vis_inter[i] = logger.vis_interior_num;
          gen_box[i] = logger.generate_box_num;
          full_box[i] = logger.full_box_num;
          skip_box[i] = logger.skip_box_num;

          // if (!std::cmp_equal(kdknn[i], query_box_seq[i].second)) {
          //   std::cout << kdknn[i] << " " << query_box_seq[i].second << " "
          //             << query_box_seq[i].first.first
          //             << query_box_seq[i].first.second << std::endl;
          //   throw std::runtime_error("wrong range query");
          // }
          assert(std::cmp_equal(kdknn[i], query_box_seq[i].second));
        });
      },
      [&]() {});

  std::cout << aveQuery << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_leaf.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(vis_inter.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(gen_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(full_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  std::cout << static_cast<double>(parlay::reduce(skip_box.cut(0, rec_num))) /
                   static_cast<double>(rec_num)
            << " " << std::flush;
  return;
}

template <typename Point, typename Tree>
typename Tree::Box ComputeBoundingBox(size_t idx, auto const& leaf_seq,
                                      auto const& inner_tree_seq) {
  using Leaf = typename Tree::Leaf;
  using Interior = typename Tree::Interior;
  using Box = typename Tree::Box;
  if (inner_tree_seq[idx].is_leaf) {
    auto const& node = inner_tree_seq[idx];
    // std::cout << node.leaf_offset << " " << node.size << std::endl;
    return Tree::GetBox(
        leaf_seq.cut(node.leaf_offset, node.leaf_offset + node.size));
  }

  Box lb, rb;
  parlay::par_do_if(
      inner_tree_seq.size() > Tree::kSerialBuildCutoff,
      [&]() {
        lb = ComputeBoundingBox<Point, Tree>(idx * 2, leaf_seq, inner_tree_seq);
      },
      [&]() {
        rb = ComputeBoundingBox<Point, Tree>(idx * 2 + 1, leaf_seq,
                                             inner_tree_seq);
      });
  return Tree::GetBox(lb, rb);
}

template <typename Point, typename Tree>
void TestBoundingBox(auto const& leaf_seq, auto const& inner_tree_seq,
                     int const& kRounds, Tree& tree) {
  using Box = typename Tree::Box;
  size_t n = leaf_seq.size();
  Box base_box, array_box, pointer_box;

  auto array_box_time = time_loop(
      kRounds, 0.1, []() {},
      [&]() {
        array_box =
            ComputeBoundingBox<Point, Tree>(1, leaf_seq, inner_tree_seq);
      },
      []() {});

  auto pointer_box_time = time_loop(
      kRounds, 0.1, []() {},
      [&]() {
        pointer_box =
            Tree::template GetBox<typename Tree::Leaf, typename Tree::Interior>(
                tree.GetRoot());
      },
      []() {});

  auto base_box_time = time_loop(
      kRounds, 0.1, []() {}, [&]() { base_box = Tree::GetBox(leaf_seq); },
      []() {});

  std::cout << "Array-tree: " << array_box_time << " " << array_box.first << " "
            << array_box.second << " " << std::endl;
  std::cout << "Pointer-tree: " << pointer_box_time << " " << pointer_box.first
            << " " << pointer_box.second << " " << std::endl;
  std::cout << "Par-for: " << base_box_time << " " << base_box.first << " "
            << base_box.second << " " << std::endl;
  return;
}

int main(int argc, char* argv[]) {
  commandLine params(
      argc, argv,
      "[-k {1,...,100}] [-d {2,3,5,7,9,10}] [-n <node num>] [-t "
      "<parallelTag>] [-p <inFile>] [-r {1,...,5}] [-q {0,1}] [-i "
      "<_insertFile>] [-s <kSummary>]");

  char* input_file_path = params.getOptionValue("-p");
  int K = params.getOptionIntValue("-k", 10);
  int dims = params.getOptionIntValue("-d", 3);
  size_t N = params.getOptionLongValue("-n", -1);
  int tag = params.getOptionIntValue("-t", 1);
  int rounds = params.getOptionIntValue("-r", 3);
  int query_type = params.getOptionIntValue("-q", 0);
  int read_insert_file = params.getOptionIntValue("-i", 1);
  int summary = params.getOptionIntValue("-s", 0);
  int tree_type = params.getOptionIntValue("-T", 0);
  int split_type = params.getOptionIntValue("-l", 0);
  int leaf_wrap = params.getOptionIntValue("-leaf", 32);

  auto test_func = []<class TreeDesc, typename Point>(
                       int const& kDim, parlay::sequence<Point> const& wp,
                       parlay::sequence<Point> const& wi, size_t const& N,
                       int const& K, int const& kRounds,
                       string const& kInsertFile, int const& kTag,
                       int const& kQueryType, int const kSummary) {
    using Tree = TreeDesc::TreeType;
    using Points = typename Tree::Points;
    using Box = typename Tree::Box;
    using InnerTreeSeq = parlay::sequence<InnerNode<Tree>>;

    Tree tree;
    constexpr bool kTestTime = true;

    BuildTree<Point, Tree, kTestTime, 2>(wp, kRounds, tree);

    Points leaf_seq(wp.size());
    size_t leaf_offset = 0;
    InnerTreeSeq inner_tree_seq(8 * wp.size());
    LinearDecomposition<Tree>(tree.GetRoot(), leaf_seq, inner_tree_seq,
                              leaf_offset, 1);
    std::cout << std::endl;
    assert(leaf_offset == wp.size());

    TestBoundingBox<Point, Tree>(leaf_seq, inner_tree_seq, kRounds, tree);
    tree.DeleteTree();
    return;

    // NOTE: knn
    Typename* kdknn;
    auto run_batch_knn = [&]<bool kIsLinear>(Points const& query_pts, int kth) {
      kdknn = new Typename[query_pts.size()];
      if constexpr (kIsLinear) {
        QueryKNNLinear<Point, Tree>(leaf_seq, inner_tree_seq, kDim, query_pts,
                                    kRounds, kdknn, kth);
      } else {
        QueryKNN<Point>(kDim, query_pts, kRounds, tree, kdknn, kth, true);
      }
      delete[] kdknn;
    };

    // NOTE: range query
    // generate boxes
    auto get_range_query_num_by_type = [](int const rec_type) {
      return rec_type == 0   ? kSmallRangeQueryNum
             : rec_type == 1 ? kMediumRangeQueryNum
                             : kLargeRangeQueryNum;
    };
    auto generate_query_box = [&](int rec_total_type, Points const& wp) {
      parlay::sequence<parlay::sequence<std::pair<typename Tree::Box, size_t>>>
          query_box_seq(rec_total_type);
      parlay::sequence<size_t> query_max_size(rec_total_type);
      for (int rec_type = 0; rec_type < rec_total_type; rec_type++) {
        auto [query_box, max_size] = gen_rectangles<Point, Tree, false, true>(
            get_range_query_num_by_type(rec_type), rec_type, wp, kDim);
        query_box_seq[rec_type] = query_box;
        query_max_size[rec_type] = max_size;
      }
      return std::make_pair(query_box_seq, query_max_size);
    };

    int const kRecType = 3;
    auto [query_box_seq, query_max_size] = generate_query_box(kRecType, wp);
    Points Out;

    // NOTE: run the range query
    auto run_range_query = [&]<bool kIsLinear>(int const rec_type) {
      int const rec_num = get_range_query_num_by_type(rec_type);
      kdknn = new Typename[rec_num];
      if constexpr (kIsLinear) {
        RangeQueryFixLinear<Point, Tree>(
            leaf_seq, inner_tree_seq, kdknn, kRounds, Out, rec_type, rec_num,
            kDim, query_box_seq[rec_type], query_max_size[rec_type]);
      } else {
        RangeQueryFix<Point, Tree>(tree, kdknn, kRounds, Out, rec_type, rec_num,
                                   kDim, query_box_seq[rec_type],
                                   query_max_size[rec_type]);
      }
      delete[] kdknn;
    };

    // NOTE: run all queries
    auto run_tests = [&]<int kIsLinear>() {
      int k[4] = {1, 10, 100, 1000};
      size_t batch_size = static_cast<size_t>(wp.size() * kBatchQueryRatio);
      for (int i = 0; i < 4; i++) {
        run_batch_knn.template operator()<kIsLinear>(wp.subseq(0, batch_size),
                                                     k[i]);
        run_batch_knn.template operator()<kIsLinear>(
            wp.subseq(wp.size() - batch_size, wp.size()), k[i]);
      }
      std::cout << std::flush;
      for (int rec_type = 0; rec_type < kRecType; rec_type++) {
        run_range_query.template operator()<kIsLinear>(rec_type);
      }
      std::cout << std::endl;
    };

    run_tests.template operator()<false>();
    run_tests.template operator()<true>();

    tree.DeleteTree();
    return;
  };

  Wrapper::ApplyTesters(tree_type, dims, split_type, params, test_func);

  return 0;
}
