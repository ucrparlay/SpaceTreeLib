// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <algorithm>
#include <iostream>

#include "parlay/parallel.h"
#include "parlay/primitives.h"

// IntV and IntE should be set depending on the size of the graphs
//  intV should have enough range to represent |V|
//  intE should have enough range to represent |E|
//  intE defaults to intV if not specified
using default_int_v_type = int;
using default_weight_type = float;

// **************************************************************
//    EDGE ARRAY REPRESENTATION
// **************************************************************

template <class intV = default_int_v_type>
struct edge {
	intV u;
	intV v;
	edge()
	{
	}
	edge(intV f, intV s) : u(f), v(s)
	{
	}
};

template <class intV = default_int_v_type>
struct edgeArray {
	parlay::sequence<edge<intV>> E;
	size_t numRows;
	size_t numCols;
	size_t nonZeros;
	edgeArray(parlay::sequence<edge<intV>> EE, size_t r, size_t c)
	    : E(std::move(EE)), numRows(r), numCols(c), nonZeros(E.size())
	{
	}
	edgeArray()
	{
	}
	edge<intV> operator[](size_t const i) const
	{
		return E[i];
	}
};

// **************************************************************
//    WEIGHED EDGE ARRAY
// **************************************************************

template <class intV = default_int_v_type, class Weight = default_weight_type>
struct wghEdge {
	intV u, v;
	Weight weight;
	wghEdge()
	{
	}
	wghEdge(intV _u, intV _v, Weight w) : u(_u), v(_v), weight(w)
	{
	}
};

template <class intV = default_int_v_type, class Weight = default_weight_type>
struct wghEdgeArray {
	using W = Weight;
	parlay::sequence<wghEdge<intV, W>> E;
	size_t n;
	size_t m;
	wghEdgeArray(parlay::sequence<wghEdge<intV, W>> E_, intV n)
	    : E(std::move(E_)), n(n), m(E.size())
	{
	}
	wghEdgeArray()
	{
	}
	wghEdge<intV> operator[](size_t const i) const
	{
		return E[i];
	}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intV = default_int_v_type>
struct vertex {
	intV const *neighbors;
	intV degree;
	vertex(intV const *N, intV const d) : neighbors(N), degree(d)
	{
	}
	vertex() : neighbors(NULL), degree(0)
	{
	}
};

template <class intV = default_int_v_type>
struct mod_vertex {
	intV *neighbors;
	intV degree;
	mod_vertex(intV *N, intV d) : neighbors(N), degree(d)
	{
	}
	mod_vertex() : neighbors(NULL), degree(0)
	{
	}
};

template <class intV = default_int_v_type, class intE = intV>
struct graph {
	using vertexId = intV;
	using edgeId = intE;
	using mvt_type = mod_vertex<intV>;
	using vt_type = vertex<intV>;
	parlay::sequence<intE> offsets;
	parlay::sequence<intV> edges;
	parlay::sequence<intV> degrees; // not always used
	size_t n;
	size_t m;
	size_t numVertices() const
	{
		return n;
	}
	size_t numEdges() const
	{
		if (degrees.size() == 0)
			return m;
		else {
			std::cout << "hello numEdges" << std::endl;
			auto dgs = parlay::delayed_seq<intE>(
				n, [&](size_t i) { return degrees[i]; });
			return parlay::reduce(dgs, parlay::addm<intE>());
		}
	}

	parlay::sequence<intE> const &get_offsets() const
	{
		return offsets;
	}

	void addDegrees()
	{
		degrees = parlay::tabulate(n, [&](size_t i) -> intV {
			return offsets[i + 1] - offsets[i];
		});
	}

	mvt_type operator[](size_t const i)
	{
		return mvt_type(edges.data() + offsets[i],
				(degrees.size() == 0)
					? offsets[i + 1] - offsets[i]
					: degrees[i]);
	}

	vt_type const operator[](size_t const i) const
	{
		return vt_type(edges.data() + offsets[i],
			       (degrees.size() == 0)
				       ? offsets[i + 1] - offsets[i]
				       : degrees[i]);
	}

	graph(parlay::sequence<intE> offsets_, parlay::sequence<intV> edges_,
	      size_t n)
	    : offsets(std::move(offsets_)), edges(std::move(edges_)), n(n),
	      m(edges.size())
	{
		if (offsets.size() != n + 1) {
			std::cout << "error in graph constructor" << std::endl;
		}
	}
};

// **************************************************************
//    WEIGHTED ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intV = default_int_v_type, class Weight = default_weight_type>
struct wghVertex {
	intV *neighbors;
	intV degree;
	Weight *nghWeights;
	wghVertex(intV *N, Weight *W, intV d)
	    : neighbors(N), nghWeights(W), degree(d)
	{
	}
};

template <class intV = default_int_v_type, class Weight = default_weight_type,
	  class intE = intV>
struct wghGraph {
	using vt_type = wghVertex<intV, Weight>;
	using W = Weight;
	parlay::sequence<intE> offsets;
	parlay::sequence<intV> edges;
	parlay::sequence<Weight> weights;
	size_t n;
	size_t m;
	size_t numVertices() const
	{
		return n;
	}
	size_t numEdges() const
	{
		return m;
	}
	// const parlay::sequence<intV>& get_offsets() const {
	//   return offsets;
	// }
	parlay::sequence<intV> get_offsets()
	{
		return offsets;
	}
	vt_type operator[](size_t const i)
	{
		return vt_type(edges.begin() + offsets[i],
			       weights.begin() + offsets[i],
			       offsets[i + 1] - offsets[i]);
	}

	wghGraph(parlay::sequence<intE> offsets_, parlay::sequence<intV> edges_,
		 parlay::sequence<Weight> weights_, size_t n)
	    : offsets(std::move(offsets_)), edges(std::move(edges_)),
	      weights(std::move(weights_)), n(n), m(edges.size())
	{
		if (offsets.size() != n + 1 || weights.size() != edges.size()) {
			std::cout << "error in weighted graph constructor"
				  << std::endl;
		}
	}
};

template <typename intV>
struct flow_graph {
	wghGraph<intV> g;
	intV source, sink;
	flow_graph(wghGraph<intV> g, intV source, intV sink)
	    : g(g), source(source), sink(sink)
	{
	}
	flow_graph copy()
	{
		return flow_graph(g.copy(), source, sink);
	}
	void del()
	{
		g.del();
	}
};
