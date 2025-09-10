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
    InnerTreeSeq inner_tree_seq(wp.size());
    size_t leaf_offset = 0;
    LinearDecomposition<Tree>(tree.GetRoot(), leaf_seq, inner_tree_seq,
                              leaf_offset, 1);
    assert(leaf_offset == wp.size());

    // begin range query test
    auto generate_query_box = [&](int rec_num, int rec_total_type,
                                  Points const& wp) {
      // NOTE: generate rectangles for the first half of the points
      parlay::sequence<parlay::sequence<std::pair<typename Tree::Box, size_t>>>
          query_box_seq(rec_total_type);
      parlay::sequence<size_t> query_max_size(rec_total_type);
      for (int i = 0; i < rec_total_type; i++) {
        auto [query_box, max_size] =
            gen_rectangles<Point, Tree, false, true>(rec_num, i, wp, kDim);
        query_box_seq[i] = query_box;
        query_max_size[i] = max_size;
      }
      return std::make_pair(query_box_seq, query_max_size);
    };

    auto [query_box_seq, query_max_size] =
        generate_query_box(kRangeQueryNum, 3, wp);
    Typename* kdknn = new Typename[kRangeQueryNum];
    Points Out;

    puts("");
    for (int rec_type = 0; rec_type < 3; rec_type++) {
      puts("-- kdtree pointer based:");
      RangeQueryFix<Point, Tree>(tree, kdknn, kRounds, Out, rec_type,
                                 kRangeQueryNum, kDim, query_box_seq[rec_type],
                                 query_max_size[rec_type]);
      std::cout << std::endl;

      puts(">> kdtree array based:");
      RangeQueryFixLinear<Point, Tree>(leaf_seq, inner_tree_seq, kdknn, kRounds,
                                       Out, rec_type, kRangeQueryNum, kDim,
                                       query_box_seq[rec_type],
                                       query_max_size[rec_type]);
      std::cout << std::endl;
    }

    delete[] kdknn;
    tree.DeleteTree();
    return;
  };

  Wrapper::ApplyTesters(tree_type, dims, split_type, params, test_func);

  return 0;
}
