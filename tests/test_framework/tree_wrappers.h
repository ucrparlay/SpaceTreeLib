#pragma once

#include <cstdint>
#include <iostream>
#include <string>

#include "../common/parse_command_line.h"
#include "augmentation_types.h"
#include "baselines/cpam_raw/cpamtree.hpp"
#include "baselines/zdtree/zdtree.hpp"
#include "dependence/concepts.h"
#include "parlay/internal/group_by.h"
#include "psi/dependence/splitter.h"
#include "psi/kd_tree.h"
#include "psi/orth_tree.h"
#include "psi/p_tree.h"
#include "test_config.h"
#include "test_functions.h"
#include "test_helpers.h"

using namespace psi;

template <typename TreeWrapper>
void PrintTreeParam() {
  std::cout << "Tree: " << TreeWrapper::TreeType::GetTreeName() << "; "
            << "AugType: " << TreeWrapper::TreeType::CheckHasBox() << "; "
            << "Split: " << TreeWrapper::SplitRule::GetSplitName() << "; "
            << "BDO: " << TreeWrapper::TreeType::GetBuildDepthOnce() << "; "
            << "Inba: " << TreeWrapper::TreeType::GetImbalanceRatio() << "; "
            << "LeafWrap: "
            << static_cast<int>(TreeWrapper::TreeType::kLeaveWrap) << "; ";

  if constexpr (std::is_integral_v<typename TreeWrapper::Point::Coord>) {
    std::cout << "Coord: integer"
              << "; ";
  } else if (std::is_floating_point_v<typename TreeWrapper::Point::Coord>) {
    std::cout << "Coord: float"
              << "; ";
  }
  std::cout << "\n" << std::flush;
  return;
}

// NOTE: default test functions for all custom tree
static auto constexpr DefaultTestFunc = []<class TreeDesc, typename Point>(
                                            int const& kDim,
                                            parlay::sequence<Point> const& wp,
                                            parlay::sequence<Point> const& wi,
                                            size_t const& N, int const& K,
                                            int const& kRounds,
                                            std::string const& kInsertFile,
                                            int const& kTag,
                                            int const& kQueryType,
                                            int const kSummary) {
  using Tree = TreeDesc::TreeType;
  using Points = typename Tree::Points;

  Tree tree(wp.size());
  constexpr bool kTestTime = true;
  BuildTree<Point, Tree, kTestTime, 0>(wp, kRounds, tree);

  // NOTE: batch insert
  if (kTag & (1 << 0)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, ratios[i]);
      }
    } else {
      BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, 0.1);
    }
  }

  // NOTE: batch delete
  if (kTag & (1 << 1)) {
    if (kSummary) {
      parlay::sequence<double> const ratios = {0.0001, 0.001, 0.01, 0.1};
      for (size_t i = 0; i < ratios.size(); i++) {
        BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds, ratios[i]);
      }
    } else {
      BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds,
                                          kBatchInsertRatio);
    }
  }

  Typename* kdknn = nullptr;
  auto run_batch_knn = [&](Points const& query_pts, int kth) {
    kdknn = new Typename[query_pts.size()];
    QueryKNN<Point>(kDim, query_pts, kRounds, tree, kdknn, kth, true);
    delete[] kdknn;
  };

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

  auto incre_update_test_bundle = [&](auto const& query_box_seq,
                                      auto const& query_max_size) {
    // NOTE: knn query
    {
      int k[3] = {1, 10, 100};

      std::cout << "in-dis-skewed knn time: ";
      size_t batch_size = static_cast<size_t>(wp.size() * kBatchQueryRatio);
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp.subseq(0, batch_size), k[i]);
      }
      puts("");

      std ::cout << "out-dis-skewed knn time: ";
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp.subseq(wp.size() - batch_size, wp.size()), k[i]);
      }
      puts("");

      // NOTE: sample points within the whole input datasets
      auto query_pts =
          parlay::pack(wp, parlay::tabulate(wp.size(), [&](size_t i) -> bool {
                         return i % (wp.size() / (batch_size * 2)) == 0;
                       }));

      std::cout << "in-dis-uniform knn time: ";
      for (int i = 0; i < 3; i++) {
        run_batch_knn(parlay::random_shuffle(query_pts.subseq(0, batch_size)),
                      k[i]);
      }
      puts("");

      std::cout << "out-dis-uniform knn time: ";
      for (int i = 0; i < 3; i++) {
        run_batch_knn(parlay::random_shuffle(
                          query_pts.subseq(batch_size, query_pts.size())),
                      k[i]);
      }
      puts("");
    }

    // NOTE: range count
    {
      int rec_num = query_box_seq[0].size();
      kdknn = new Typename[rec_num];

      std::cout << "range count time: ";
      for (int i = 0; i < 3; i++) {
        RangeCountFix<Point>(tree, kdknn, kRounds, i, rec_num, kDim,
                             query_box_seq[i], query_max_size[i]);
      }
      delete[] kdknn;
      puts("");
    }

    // NOTE: range query
    {
      int rec_num = query_box_seq[0].size();
      kdknn = new Typename[rec_num];

      std::cout << "range query time: ";
      for (int i = 0; i < 3; i++) {
        Points Out;
        RangeQueryFix<Point>(tree, kdknn, kRounds, Out, i, rec_num, kDim,
                             query_box_seq[i], query_max_size[i]);
      }
      delete[] kdknn;
      puts("");
    }
  };

  // NOTE: scalability
  if (kTag & (1 << 2)) {
    puts("");
    BuildTree<Point, Tree, kTestTime, 0>(wp, kRounds, tree);
    BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds,
                                        kBatchInsertRatio);
    BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds,
                                        kBatchInsertRatio);
  }

  // NOTE: batch insert by step
  if (kTag & (1 << 3)) {
    puts("");
    BuildTree<Point, Tree, kTestTime, 3>(wp, kRounds, tree, 2);

    auto [query_box_seq, query_max_size] =
        generate_query_box(kRangeQueryNum, 3, wp.subseq(0, wp.size() / 2));

    incre_update_test_bundle(query_box_seq, query_max_size);

    parlay::sequence<double> const ratios = {0.1, 0.01, 0.001, 0.0001};
    for (auto rat : ratios) {
      BatchInsertByStep<Point, Tree, true>(tree, wp, kRounds, rat);
      incre_update_test_bundle(query_box_seq, query_max_size);
    }
  }

  // NOTE: batch delete by step
  if (kTag & (1 << 4)) {
    puts("");
    auto [query_box_seq, query_max_size] =
        generate_query_box(kRangeQueryNum, 3, wp.subseq(0, wp.size() / 2));
    BuildTree<Point, Tree, kTestTime, 3>(wp, kRounds, tree, 2);
    incre_update_test_bundle(query_box_seq, query_max_size);

    parlay::sequence<double> const ratios = {0.1, 0.01, 0.001, 0.0001};
    for (auto rat : ratios) {
      BatchDeleteByStep<Point, Tree, true>(tree, wp, kRounds, rat);
      incre_update_test_bundle(query_box_seq, query_max_size);
    }
  }

  // real world
  if (kTag & (1 << 5)) {
    puts("");
    auto [query_box_seq, query_max_size] =
        generate_query_box(kRangeQueryNum, 3, wp.subseq(0, wp.size() / 2));

    BuildTree<Point, Tree, kTestTime, 3>(wp, kRounds, tree, 2);
    incre_update_test_bundle(query_box_seq, query_max_size);

    BatchInsertByStep<Point, Tree, true>(tree, wp, kRounds, 0.0001);
    incre_update_test_bundle(query_box_seq, query_max_size);

    BatchDeleteByStep<Point, Tree, true>(tree, wp, kRounds, 0.0001);
    incre_update_test_bundle(query_box_seq, query_max_size);
  }

  // range query with log
  if (kTag & (1 << 6)) {
    puts("");
    BatchInsertByStep<Point, Tree, true>(tree, wp, kRounds, 0.0001);

    auto [query_box_seq, query_max_size] = generate_query_box(
        kSingleQueryLogRepeatNum, 3, wp.subseq(0, wp.size() / 2));

    // NOTE: range query
    {
      int rec_num = kSingleQueryLogRepeatNum;
      kdknn = new Typename[rec_num];

      std::cout << "range query time: " << std::endl;
      for (int i = 0; i < 3; i++) {
        Points Out;
        RangeQuerySerialWithLog<Point>(tree, kdknn, kRounds, Out, i, rec_num,
                                       kDim, query_box_seq[i],
                                       query_max_size[i]);
      }
      delete[] kdknn;
      puts("");
    }
  }

  // range query with log
  if (kTag & (1 << 7)) {
    puts("");
    parlay::sequence<double> const ratios = {
        0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01,
        0.02,   0.05,   0.1,    0.2,   0.5,   1.0};
    std::cout << "Insert: ";
    for (size_t i = 0; i < ratios.size(); i++) {
      BatchInsert<Point, Tree, kTestTime>(tree, wp, wi, kRounds, ratios[i]);
    }
    std::cout << std::endl;
    std::cout << "Delete: ";
    for (size_t i = 0; i < ratios.size(); i++) {
      BatchDelete<Point, Tree, kTestTime>(tree, wp, wp, kRounds, ratios[i]);
    }
    puts("");
  }

  if (kQueryType & (1 << 0)) {  // NOTE: KNN
    size_t batch_size = static_cast<size_t>(wp.size() * kBatchQueryRatio);

    if (kSummary == 0) {
      int k[3] = {1, 10, 100};
      for (int i = 0; i < 3; i++) {
        run_batch_knn(wp.subseq(0, batch_size), k[i]);
      }
    } else {  // test kSummary
      run_batch_knn(wp.subseq(0, batch_size), K);
    }
  }

  if (kQueryType & (1 << 1)) {  // NOTE: range count
    auto [query_box_seq, query_max_size] =
        generate_query_box(kRangeQueryNum, 3, wp);

    if (!kSummary) {
      int recNum = kRangeQueryNum;
      kdknn = new Typename[recNum];

      for (int i = 0; i < 3; i++) {
        RangeCountFix<Point>(tree, kdknn, kRounds, i, recNum, kDim,
                             query_box_seq[i], query_max_size[i]);
      }

      delete[] kdknn;
    }
  }

  if (kQueryType & (1 << 2)) {  // NOTE: range query
    if (kSummary == 0) {
      auto [query_box_seq, query_max_size] =
          generate_query_box(kRangeQueryNum, 3, wp);

      int recNum = kRangeQueryNum;
      kdknn = new Typename[recNum];

      for (int i = 0; i < 3; i++) {
        Points Out;
        RangeQueryFix<Point>(tree, kdknn, kRounds, Out, i, recNum, kDim,
                             query_box_seq[i], query_max_size[i]);
      }
      delete[] kdknn;
    } else if (kSummary == 1) {  // NOTE: for kSummary
      auto [query_box_seq, query_max_size] =
          generate_query_box(kSummaryRangeQueryNum, 3, wp);

      kdknn = new Typename[kSummaryRangeQueryNum];
      Points Out;
      RangeQueryFix<Point>(tree, kdknn, kRounds, Out, 2, kSummaryRangeQueryNum,
                           kDim, query_box_seq[2], query_max_size[2]);
      delete[] kdknn;
    }
  }

  std::cout << "\n" << std::flush;

  tree.DeleteTree();

  return;
};

class Wrapper {
 public:
  // NOTE: determine the build depth once for the orth tree
  static consteval uint8_t OrthGetBuildDepthOnce(uint8_t const dim) {
    if (dim == 2 || dim == 3) {
      return 6;
    } else if (dim == 4) {
      return 8;
    } else if (dim >= 5 && dim <= 8) {
      return dim;
    } else {
      static_assert("Cannot decide the build tree depth once for this dim");
      return 0;
    }
  }

  // NOTE: Trees
  template <class PointType, class SplitRuleType, class LeafAugType,
            class InteriorAugType>
  struct KdTreeWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename psi::KdTree<
        psi::TypeTrait<Point, psi::pointer_view::Node, SplitRule, LeafAugType,
                       InteriorAugType>>;
  };

  template <class PointType, class SplitRuleType, class InteriorAugType>
  struct ArrayTreeWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename psi::array_view::KdTreeArray<psi::TypeTrait<
        Point, uint32_t, SplitRule, InteriorAugType, InteriorAugType, 1>>;
    // BUG: using uint32_t as the node type
  };

  template <class PointType, class SplitRuleType, class LeafAugType,
            class InteriorAugType>
  struct OrthTreeWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename psi::OrthTree<psi::TypeTrait<
        Point, psi::pointer_view::Node, SplitRule, LeafAugType, InteriorAugType,
        Point::GetDim(), OrthGetBuildDepthOnce(Point::GetDim())>>;
  };

  template <class PointType, class SplitRuleType>
  struct PTreeWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename psi::PTree<
        psi::TypeTrait<Point, psi::pointer_view::Node, SplitRule>>;
  };

  template <class PointType, class SplitRuleType>
  struct CpamRawWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename CPAMTree::CpamRaw<
        psi::TypeTrait<Point, psi::pointer_view::Node, SplitRule>>;
  };

  /* Zdtree Wrapper */
  template <class PointType, class SplitRuleType>
  struct ZdTreeWrapper {
    using Point = PointType;
    using SplitRule = SplitRuleType;
    using TreeType = typename ZD::Zdtree<
        psi::TypeTrait<Point, psi::pointer_view::Node, SplitRule>>;
  };

  // NOTE: Apply the dim and split rule
  struct AugId {
    using IdType = int;
    IdType id;

    bool operator<(AugId const& rhs) const { return id < rhs.id; }
    bool operator==(AugId const& rhs) const { return id == rhs.id; }
    friend std::ostream& operator<<(std::ostream& os, AugId const& rhs) {
      os << rhs.id;
      return os;
    }
  };

  // For the spacial filling curve, we use the AugIdCode to
  // ensure the id is unique and the code is used to determine the
  // order of the points in the tree.
  struct AugIdCode {
    using IdType = int_fast32_t;
    using CurveCode = uint64_t;

    AugIdCode() : code(0), id(0) {}

    void SetMember(CurveCode const& val) { code = val; }

    bool operator<(AugIdCode const& rhs) const {
      return code == rhs.code ? id < rhs.id : code < rhs.code;
    }

    bool operator==(AugIdCode const& rhs) const {
      // return code == rhs.code && id == rhs.id;
      // WARN: code is not important, we only need to ensure the id
      return id == rhs.id;
    }

    friend std::ostream& operator<<(std::ostream& os, AugIdCode const& rhs) {
      os << rhs.code << " " << rhs.id;
      return os;
    }

    CurveCode code;
    IdType id;
  };

  // NOTE: driven functions
  template <typename TreeWrapper, typename RunFunc>
  static void Run(commandLine& P, RunFunc test_func) {
    char* input_file_path = P.getOptionValue("-p");
    int K = P.getOptionIntValue("-k", 100);
    int dims = P.getOptionIntValue("-d", 3);
    size_t N = P.getOptionLongValue("-n", -1);
    int tag = P.getOptionIntValue("-t", 1);
    int rounds = P.getOptionIntValue("-r", 3);
    int query_type = P.getOptionIntValue("-q", 0);
    int read_insert_file = P.getOptionIntValue("-i", 1);
    char* insert_file_path_cml = P.getOptionValue("-I");
    int summary = P.getOptionIntValue("-s", 0);
    int tree_type = P.getOptionIntValue("-T", 0);
    int split_type = P.getOptionIntValue("-l", 0);

    using Point = typename TreeWrapper::Point;
    using Points = parlay::sequence<Point>;
    constexpr auto kDim = Point::GetDim();

    PrintTreeParam<TreeWrapper>();

    std::string name, insert_file_path = "";
    Points wp, wi;

    if (input_file_path != NULL) {  // NOTE: read main Points
      name = std::string(input_file_path);
      name = name.substr(name.rfind('/') + 1);
      std::cout << name << " ";
      auto [n, d] = read_points<Point>(input_file_path, wp, 0);
      N = n;
      assert(d == kDim);
    }

    if (read_insert_file == 1) {           // NOTE: read points to be inserted
      if (insert_file_path_cml != NULL) {  // given in commadnline
        insert_file_path = std::string(insert_file_path_cml);
      } else {  // determine the name otherwise
        int id = std::stoi(name.substr(0, name.find_first_of('.')));
#ifdef CCP
        id = (id + 1) % 10;  // WARN: MOD graph number used to test
#else
        id = (id + 1) % 3;
#endif  // CCP
        if (!id) id++;
        auto pos = std::string(input_file_path).rfind('/') + 1;
        insert_file_path = std::string(input_file_path).substr(0, pos) +
                           std::to_string(id) + ".in";
      }
      auto [n, d] = read_points<Point>(insert_file_path.c_str(), wi, N);
      assert(d == kDim);
    }

    // Apply the test function
    test_func.template operator()<TreeWrapper, Point>(
        kDim, wp, wi, N, K, rounds, insert_file_path, tag, query_type, summary);
  };

  // NOTE: For kd tree and orth tree
  template <typename RunFunc>
  static void ApplyOrthogonal(int const tree_type, int const dim,
                              int const split_type, commandLine& params,
                              RunFunc test_func) {
    auto build_tree_type = [&]<typename Point, typename SplitRule>() {
      using BT = psi::BaseTree<psi::TypeTrait<Point>>;
      if (tree_type == 0) {
        Run<KdTreeWrapper<Point, SplitRule, LeafAugBox<BT>,
                          InteriorAugBox<BT>>>(params, test_func);
      } else if (tree_type == 1) {
        Run<OrthTreeWrapper<Point, SplitRule, LeafAugBox<BT>,
                            InteriorAugBox<BT>>>(params, test_func);
      } else if (tree_type == 4) {  // NOTE: for boost
        Run<KdTreeWrapper<Point, SplitRule, LeafAugBox<BT>,
                          InteriorAugBox<BT>>>(params, test_func);
      }
    };

    // NOTE: pick the split rule
    // The lsb is the dim rule and the msb is the divide rule
    auto run_with_split_type = [&]<typename Point>() {
      if (!(split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 0 -> max_stretch + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 1 -> rotate_dim + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if (!(split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 2 -> max_stretch + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 3 -> rotate + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else {
        std::cout << "Unsupported split type: " << split_type << std::endl;
      }
    };

    if (dim == 2) {
      run_with_split_type.template operator()<AugPoint<Coord, 2, AugId>>();
    }
  }

  template <typename RunFunc>
  static void ApplySpacialFillingCurve(int const tree_type, int const dim,
                                       int const split_type,
                                       commandLine& params, RunFunc test_func) {
    auto build_tree_type = [&]<typename Point, typename SplitRule>() {
      if (tree_type == 0) {
        // run.template operator()<KdTreeWrapper<Point, SplitRule>>();
      } else if (tree_type == 1) {
        // run.template operator()<OrthTreeWrapper<Point, SplitRule>>();
      } else if (tree_type == 2) {
        Run<PTreeWrapper<Point, SplitRule>>(params, test_func);
      } else {
        std::cout << "Unsupported tree type: " << tree_type << std::endl;
      }
    };

    // NOTE: pick the split rule
    auto run_with_split_type = [&]<typename Point>() {
      if (split_type & (1 << 0)) {
        build_tree_type.template operator()<
            Point,
            psi::SpacialFillingCurve<HilbertCurve<psi::TypeTrait<Point>>>>();
      } else if (split_type & (1 << 1)) {
        build_tree_type.template operator()<
            Point,
            psi::SpacialFillingCurve<MortonCurve<psi::TypeTrait<Point>>>>();
      }
    };

    if (dim == 2) {
      run_with_split_type.template operator()<AugPoint<Coord, 2, AugIdCode>>();
    }
  }

  template <typename RunFunc>
  static void ApplyBaselines(int const tree_type, int const dim,
                             int const split_type, commandLine& params,
                             RunFunc test_func) {
    auto build_tree_type = [&]<typename Point, typename SplitRule>() {
      if (tree_type == 0) {
        // run.template operator()<KdTreeWrapper<Point, SplitRule>>();
      } else if (tree_type == 1) {
        // run.template operator()<OrthTreeWrapper<Point, SplitRule>>();
      } else if (tree_type == 2) {
        // Run<PTreeWrapper<Point, SplitRule>>(params, test_func);
      } else if (tree_type == 3) {
        Run<CpamRawWrapper<Point, SplitRule>>(params, test_func);
      } else if (tree_type == 4) {
        ;  // for boost
      } else if (tree_type == 5) {
        Run<ZdTreeWrapper<typename ZD::geobase::Point, SplitRule>>(params,
                                                                   test_func);
      } else {
        std::cout << "Unsupported tree type: " << tree_type << std::endl;
      }
    };

    // NOTE: pick the split rule
    auto run_with_split_type = [&]<typename Point>() {
      if (split_type & (1 << 0)) {
        build_tree_type.template operator()<
            Point,
            psi::SpacialFillingCurve<HilbertCurve<psi::TypeTrait<Point>>>>();
      } else if (split_type & (1 << 1)) {
        build_tree_type.template operator()<
            Point,
            psi::SpacialFillingCurve<MortonCurve<psi::TypeTrait<Point>>>>();
      }
    };

    if (dim == 2) {
      run_with_split_type.template operator()<AugPoint<Coord, 2, AugIdCode>>();
    } else if (dim == 3) {
      if (tree_type == 4) {
        ;
      } else {
        run_with_split_type
            .template operator()<AugPoint<Coord, 3, AugIdCode>>();
      }
    }
  }

  template <typename RunFunc>
  static void ApplyTesters(int const tree_type, int const dim,
                           int const split_type, commandLine& params,
                           RunFunc test_func) {
    auto build_tree_type = [&]<typename Point, typename SplitRule>() {
      using BT = psi::BaseTree<psi::TypeTrait<Point>>;
      if (tree_type == 0) {
        Run<KdTreeWrapper<Point, SplitRule, LeafAugBox<BT>,
                          InteriorTester<BT>>>(params, test_func);
      }
    };

    // NOTE: pick the split rule
    auto run_with_split_type = [&]<typename Point>() {
      if (!(split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 0 -> max_stretch + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 1 -> rotate_dim + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if (!(split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 2 -> max_stretch + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 3 -> rotate + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else {
        std::cout << "Unsupported split type: " << split_type << std::endl;
      }
    };

    if (dim == 2) {
      run_with_split_type.template operator()<AugPoint<Coord, 2, AugIdCode>>();
    } else if (dim == 3) {
      run_with_split_type.template operator()<AugPoint<Coord, 3, AugIdCode>>();
    }
  }

  template <typename RunFunc>
  static void ApplyArrayTree(int const tree_type, int const dim,
                             int const split_type, commandLine& params,
                             RunFunc test_func) {
    auto build_tree_type = [&]<typename Point, typename SplitRule>() {
      using BT = psi::BaseTree<psi::TypeTrait<Point>>;
      if (tree_type == 0) {
        Run<ArrayTreeWrapper<Point, SplitRule, NodeAugBox<BT>>>(params,
                                                                test_func);
      }
    };

    // NOTE: pick the split rule
    // The lsb is the dim rule and the msb is the divide rule
    auto run_with_split_type = [&]<typename Point>() {
      if (!(split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 0 -> max_stretch + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && !(split_type & (1 << 1))) {
        // NOTE: 1 -> rotate_dim + object_mid
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::ObjectMedian<psi::TypeTrait<Point>>>>();
      } else if (!(split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 2 -> max_stretch + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::MaxStretchDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else if ((split_type & (1 << 0)) && (split_type & (1 << 1))) {
        // NOTE: 3 -> rotate + spatial_median
        build_tree_type.template
        operator()<Point, psi::OrthogonalSplitRule<
                              psi::RotateDim<psi::TypeTrait<Point>>,
                              psi::SpatialMedian<psi::TypeTrait<Point>>>>();
      } else {
        std::cout << "Unsupported split type: " << split_type << std::endl;
      }
    };

    if (dim == 2) {
      run_with_split_type.template operator()<AugPoint<Coord, 2, AugId>>();
    }
  }
};
