#pragma once

#include <parlay/parallel.h>
#include "../base_tree.h"

namespace cpdd {

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
bool BaseTree<Point, kBDO>::CheckBox(Node* T, const Box& bx) {
    assert(T != nullptr);
    assert(LegalBox(bx));
    Points wx = Points::uninitialized(T->size);
    LOG << T->size << ENDL;
    FlattenRec<Leaf, Interior>(T, parlay::make_slice(wx));
    auto b = GetBox(parlay::make_slice(wx));
    // LOG << b.first << b.second << ENDL;
    return WithinBox(b, bx);
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
size_t BaseTree<Point, kBDO>::CheckSize(Node* T) {
    if (T->is_leaf) {
        return T->size;
    }
    if constexpr (IsBinaryNode<Interior>) {
        Interior* TI = static_cast<Interior*>(T);
        size_t l = CheckSize<Leaf, Interior>(TI->left);
        size_t r = CheckSize<Leaf, Interior>(TI->right);
        assert(l + r == T->size);
        return T->size;
    } else {
        assert(IsMultiNode<Interior>);
        Interior* TI = static_cast<Interior*>(T);
        size_t sum = 0;
        for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
            sum += CheckSize<Leaf, Interior>(TI->tree_nodes[i]);
        }
        assert(sum == T->size);
        return T->size;
    }
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
void BaseTree<Point, kBDO>::CheckTreeSameSequential(Node* T, int dim,
                                                    const int& DIM) {
    if (T->is_leaf) {
        // assert( PickRebuildDim( T, DIM ) == dim );
        return;
    }
    if constexpr (IsBinaryNode<Interior>) {
        Interior* TI = static_cast<Interior*>(T);
        if (TI->split.second != dim) {
            LOG << int(TI->split.second) << " " << int(dim) << " " << TI->size
                << ENDL;
        }
        assert(TI->split.second == dim);
        dim = (dim + 1) % DIM;
        parlay::par_do_if(
            T->size > 1000,
            [&]() {
                CheckTreeSameSequential<Leaf, Interior>(TI->left, dim, DIM);
            },
            [&]() {
                CheckTreeSameSequential<Leaf, Interior>(TI->right, dim, DIM);
            });
    } else {
        assert(IsMultiNode<Interior>);
        Interior* TI = static_cast<Interior*>(T);
        assert((1 << TI->split.size()) == TI->tree_nodes.size());
        for (int i = 0; i < TI->split.size(); i++) {
            assert(TI->split[i].second == dim);
            dim += 1;
        }
        assert(dim == DIM);
        for (int i = 0; i < TI->tree_nodes.size(); i++) {
            CheckTreeSameSequential<Leaf, Interior>(TI->tree_nodes[i], 0, DIM);
        }
    }
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior, typename SplitRule>
void BaseTree<Point, kBDO>::Validate(const DimsType DIM) {
    if (CheckBox<Leaf, Interior>(this->root_, this->tree_box_) &&
        LegalBox(this->tree_box_)) {
        std::cout << "Correct bounding Box" << std::endl << std::flush;
    } else {
        std::cout << "wrong bounding Box" << std::endl << std::flush;
        abort();
    }

    // NOTE: used to check rotate dimension
    // For kdtree binary node, the dummy node may break the rotation manner
    if constexpr (IsRotateDimSplit<SplitRule> && IsMultiNode<Interior>) {
        CheckTreeSameSequential<Leaf, Interior>(this->root_, 0, DIM);
        std::cout << "Correct rotate dimension" << std::endl << std::flush;
    }

    if (CheckSize<Leaf, Interior>(this->root_) == this->root_->size) {
        std::cout << "Correct size" << std::endl << std::flush;
    } else {
        std::cout << "wrong tree size" << std::endl << std::flush;
        abort();
    }
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
size_t BaseTree<Point, kBDO>::GetTreeHeight() {
    size_t deep = 0;
    return GetMaxTreeDepth<Leaf, Interior>(this->root_, deep);
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
size_t BaseTree<Point, kBDO>::GetMaxTreeDepth(Node* T, size_t deep) {
    if (T->is_leaf) {
        return deep;
    }

    Interior* TI = static_cast<Interior*>(T);
    if constexpr (IsBinaryNode<Interior>) {
        int l = GetMaxTreeDepth<Leaf, Interior>(TI->left, deep + 1);
        int r = GetMaxTreeDepth<Leaf, Interior>(TI->right, deep + 1);
        return std::max(l, r);
    } else {
        size_t max_depth = 0;
        for (int i = 0; i < TI->tree_nodes.size(); i++) {
            max_depth = std::max(max_depth, GetMaxTreeDepth<Leaf, Interior>(
                                                TI->tree_nodes[i], deep + 1));
        }
        return max_depth;
    }
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
double BaseTree<Point, kBDO>::GetAveTreeHeight() {
    // LOG << IsBinaryNode<Interior> << ENDL;
    parlay::sequence<size_t> heights(this->root_->size);
    size_t idx = 0;
    CountTreeHeights<Leaf, Interior>(this->root_, 0, idx, heights);
    // auto kv = parlay::histogram_by_key( heights.cut( 0, idx ) );
    // std::sort( kv.begin(), kv.end(),
    //            [&]( auto a, auto b ) { return a.first < b.first; } );
    // for ( auto i : kv )
    //     LOG << i.first << " " << i.second << ENDL;
    return double(1.0 * parlay::reduce(heights.cut(0, idx)) / idx);
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
size_t BaseTree<Point, kBDO>::CountTreeNodesNum(Node* T) {
    if (T->is_leaf) {
        return 1;
    }

    Interior* TI = static_cast<Interior*>(T);
    if constexpr (IsBinaryNode<Interior>) {
        size_t l, r;
        parlay::par_do(
            [&]() { l = CountTreeNodesNum<Leaf, Interior>(TI->left); },
            [&]() { r = CountTreeNodesNum<Leaf, Interior>(TI->right); });
        return l + r + 1;
    } else {
        size_t sum = 0;
        for (int i = 0; i < TI->tree_nodes.size(); i++) {
            sum += CountTreeNodesNum<Leaf, Interior>(TI->tree_nodes[i]);
        }
        return sum + 1;
    }
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
void BaseTree<Point, kBDO>::CountTreeHeights(
    Node* T, size_t deep, size_t& idx, parlay::sequence<size_t>& heights) {
    if (T->is_leaf) {
        heights[idx++] = deep;
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    if constexpr (IsBinaryNode<Interior>) {
        CountTreeHeights<Leaf, Interior>(TI->left, deep + 1, idx, heights);
        CountTreeHeights<Leaf, Interior>(TI->right, deep + 1, idx, heights);
    } else {
        for (int i = 0; i < TI->tree_nodes.size(); i++) {
            CountTreeHeights<Leaf, Interior>(TI->tree_nodes[i], deep + 1, idx,
                                             heights);
        }
    }
    return;
}

}  // namespace cpdd
