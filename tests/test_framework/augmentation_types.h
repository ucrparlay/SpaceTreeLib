#pragma once

#include <optional>

#include "psi/base_tree.h"

using psi::pointer_based::Node;

// Augmentation types for leaves and interior nodes

template <class BaseTree>
struct LeafAugEmpty {
  using BT = BaseTree;
  using Box = BT::Box;
  using Slice = BT::Slice;
  using Geo = BT::Geo;

  LeafAugEmpty() {};
  LeafAugEmpty(Box const& _box) {};
  LeafAugEmpty(Slice In) {};
  void UpdateAug(Slice In) { return; }
  void Reset() { return; }
};

template <class BaseTree>
struct LeafAugBox {
  using BT = BaseTree;
  using Box = BT::Box;
  using Slice = BT::Slice;
  using Geo = BT::Geo;

  LeafAugBox() : box(Geo::GetEmptyBox()) {};
  LeafAugBox(Box const& _box) : box(_box) {};
  LeafAugBox(Slice In) : box(Geo::GetBox(In)) {};

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  void UpdateAug(Slice In) {
    this->box = Geo::GetBox(In);
    return;
  }

  void Reset() {
    this->box = Geo::GetEmptyBox();
    return;
  }

  Box box;
};

template <class BaseTree>
struct InteriorAugEmpty {
  using BT = BaseTree;
  using Geo = BT::Geo;

  InteriorAugEmpty() { force_par_indicator.reset(); }
  InteriorAugEmpty(bool) { force_par_indicator.reset(); }

  // use a bool to reload default constructor
  template <typename Leaf, typename Interior>
  static bool Create(Node* l, Node* r) {
    return true;
  }

  template <typename TreeNodes>
  static bool Create(TreeNodes const& /*nodes*/) {
    return true;
  }

  void SetParallelFlag(bool const flag) {
    this->force_par_indicator.emplace(flag);
  }

  void ResetParallelFlag() { this->force_par_indicator.reset(); }

  bool GetParallelFlagIniStatus() {
    return this->force_par_indicator.has_value();
  }

  bool ForceParallel(size_t sz) const {
    return this->force_par_indicator.has_value()
               ? this->force_par_indicator.value()
               : sz > BT::kSerialBuildCutoff;
  }

  template <typename Leaf, typename Interior>
  void Update(Node*, Node*) {
    return;
  }

  template <typename TreeNodes>
  void Update(TreeNodes const& /*nodes*/) {
    return;
  }

  void Reset() { this->force_par_indicator.reset(); }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  std::optional<bool> force_par_indicator;
};

template <class BaseTree>
struct InteriorAugBox : public InteriorAugEmpty<BaseTree> {
  using BT = BaseTree;
  using Box = BT::Box;
  using Slice = BT::Slice;
  using BaseAug = InteriorAugEmpty<BT>;

  using Geo = BT::Geo;

  InteriorAugBox() : BaseAug(), box(Geo::GetEmptyBox()) {}
  InteriorAugBox(Box const& _box) : BaseAug(), box(_box) {}

  // binary create
  template <typename Leaf, typename Interior>
  static Box Create(Node* l, Node* r) {
    return Geo::GetBox(BT::template RetriveBox<Leaf, Interior>(l),
                       BT::template RetriveBox<Leaf, Interior>(r));
  }

  // multi create
  template <typename Leaf, typename Interior, typename TreeNodes>
  static Box Create(TreeNodes const& nodes) {
    Box box = Geo::GetEmptyBox();
    for (auto t : nodes) {
      box = Geo::GetBox(box, BT::template RetriveBox<Leaf, Interior>(t));
    }
    return box;
  }

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  // binary update
  template <typename Leaf, typename Interior>
  void Update(Node* l, Node* r) {
    this->box = this->Create<Leaf, Interior>(l, r);
    return;
  }

  // multi update
  template <typename Leaf, typename Interior, typename TreeNodes>
  void Update(TreeNodes const& nodes) {
    this->box = this->Create<Leaf, Interior>(nodes);
    return;
  }

  void Reset() {
    BaseAug::Reset();
    this->force_par_indicator.reset();
  }

  Box box;
};

// For testing purposes
template <typename BaseTree>
struct InteriorTester : InteriorAugBox<BaseTree> {
  using BT = BaseTree;
  using Box = BT::Box;
  using Slice = BT::Slice;

  InteriorTester() : InteriorAugBox<BaseTree>(), leaf_offset(0) {}
  InteriorTester(Box const& _box)
      : InteriorAugBox<BaseTree>(_box), leaf_offset(0) {}
  size_t leaf_offset;
};
