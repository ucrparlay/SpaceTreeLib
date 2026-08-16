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

#include <algorithm>
#include <iostream>

#include "atomics.h"
#include "geometry.h"
#include "get_time.h"
#include "parlay/hash_table.h"
#include "parlay/primitives.h"
#include "topology.h"

using parlay::hash64;
using parlay::hashtable;
using parlay::parallel_for;
using parlay::sequence;
using parlay::tabulate;

using std::cout;
using std::endl;
using std::less;
using std::pair;

using triang_t = triangle<point>;
using vertex_t = vertex<point>;
using simplex_t = simplex<point>;
using index_t = int;
using index_pair = pair<index_t, index_t>;
using edge = pair<index_pair, triang_t *>;

// Hash table to store skinny triangles
struct hashEdges {
	using k_type = index_pair;
	using eType = edge *;
	eType empty()
	{
		return NULL;
	}
	k_type getKey(eType v)
	{
		return v->first;
	}
	size_t hash(k_type s)
	{
		return hash64(s.first) + 3 * (hash64(s.second));
	}
	int cmp(k_type s1, k_type s2)
	{
		return ((s1.first > s2.first)	  ? 1
			: (s1.first < s2.first)	  ? -1
			: (s1.second > s2.second) ? 1
			: (s1.second < s2.second) ? -1
						  : 0);
	}
	bool cas(eType *p, eType o, eType n)
	{
		return pbbs::atomic_compare_and_swap(p, o, n);
	}
	bool replaceQ(eType s, eType s2)
	{
		return 0;
	}
};

using edge_table_type = hashtable<hashEdges>;

edge_table_type makeEdgeTable(size_t m)
{
	return edge_table_type(m, hashEdges());
}

std::pair<sequence<triang_t>, sequence<vertex_t>>
topology_from_triangles(triangles<point> &tris, size_t extra_points = 0)
{
	size_t n = tris.numPoints();
	size_t m = tris.numTriangles();

	auto V = tabulate(n + extra_points, [&](size_t i) {
		return (i < n) ? vertex_t(tris.P[i], i) : vertex_t();
	});

	sequence<triang_t> triangs(m + 2 * extra_points);
	sequence<edge> E(m * 3);
	edge_table_type et_type = makeEdgeTable(m * 6);
	parallel_for(0, m, [&](size_t i) {
		for (int j = 0; j < 3; j++) {
			E[i * 3 + j] = edge(index_pair(tris.T[i][j],
						       tris.T[i][(j + 1) % 3]),
					    &triangs[i]);
			et_type.insert(&E[i * 3 + j]);
			triangs[i].vtx[(j + 2) % 3] = &V[tris.T[i][j]];
		}
	});

	parallel_for(0, m, [&](size_t i) {
		triangs[i].id = i;
		triangs[i].initialized = 1;
		triangs[i].bad = 0;
		for (int j = 0; j < 3; j++) {
			index_pair key = {tris.T[i][(j + 1) % 3], tris.T[i][j]};
			edge *ed = et_type.find(key);
			if (ed != NULL)
				triangs[i].ngh[j] = ed->second;
			else {
				triangs[i].ngh[j] = NULL;
				// triangs[i].vtx[j]->boundary = 1;
				// triangs[i].vtx[(j+2)%3]->boundary = 1;
			}
		}
	});
	return std::pair(std::move(triangs), std::move(V));
}

// Note that this is not currently a complete test of correctness
// For example it would allow a set of disconnected triangles, or even no
// triangles
bool check_delaunay(sequence<triang_t> &tri_seq, size_t boundary_size)
{
	size_t n = tri_seq.size();
	sequence<size_t> boundary_count(n, 0);
	size_t insideOutError = n;
	size_t inCircleError = n;
	parallel_for(0, n, [&](size_t i) {
		if (tri_seq[i].initialized) {
			simplex_t t = simplex(&tri_seq[i], 0);
			for (int j = 0; j < 3; j++) {
				simplex_t a = t.across();
				if (a.valid()) {
					vertex_t *v =
						a.rotClockwise().firstVertex();

					// Check that the neighbor is outside
					// the triangle
					if (!t.outside(v)) {
						double vz = triAreaNormalized(
							t.t->vtx[(t.o + 2) % 3]
								->pt,
							v->pt,
							t.t->vtx[t.o]->pt);
						// allow for small error
						if (vz < -1e-10)
							pbbs::write_min(
								&insideOutError,
								i,
								less<size_t>());
					}

					// Check that the neighbor is not in
					// circumcircle of the triangle
					if (t.inCirc(v)) {
						double vz = inCircleNormalized(
							t.t->vtx[0]->pt,
							t.t->vtx[1]->pt,
							t.t->vtx[2]->pt, v->pt);
						// allow for small error
						if (vz > 1e-10)
							pbbs::write_min(
								&inCircleError,
								i,
								less<size_t>());
					}
				} else
					boundary_count[i]++;
				t = t.rotClockwise();
			}
		}
	});
	// if (boundary_size != reduce(boundary_count))
	//   cout << "Wrong boundary size: should be " << boundary_size
	// 	 << " is " << reduce(boundary_count) << endl;

	if (insideOutError < n) {
		cout << "delaunayCheck: neighbor inside triangle at triangle "
		     << inCircleError << endl;
		return 1;
	}
	if (inCircleError < n) {
		cout << "In Circle Violation at triangle " << inCircleError
		     << endl;
		return 1;
	}

	return 0;
}
