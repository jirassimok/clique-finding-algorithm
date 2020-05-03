/**
 * @file graph.h
 * @author Jacob Komissar
 * @date 2020-05-02
 *
 * Graph macro definitions.
 *
 * This file declares a constant-sized graph and macros that operate on it.
 */
#ifndef GRAPH_H
#define GRAPH_H

#include <stdbool.h>

/**
 * Number of vertices in the graph.
 *
 * The program may fail for N < 2.
 */
#define N (1024)

/**
 * Subgraphs as logical arrays: subgraph[u] indicates whether vertex u is in
 * the subgraph.
 */
typedef bool subgraph[N];

typedef int vertex;

/**
 * The adjacency matrix
 */
static bool ADJ[N][N];

/**
 * The whole graph represented as a subgraph, suitable for use with the
 * subgraph macros.
 */
static const subgraph GRAPH = {[0 ... N-1] = true};

/**
 * Declare a macro that acts differently with 2 or 3 arguments.
   #define MACRO(...) SHIFT_MACRO(__VA_ARGS__, m2args, m1arg)(__VA_ARGS__)
 * To use for 1 or 2 arguments, pass any value as an additional first argument.
 */
#define SHIFT_MACRO(_1, _2, _3, ...) SHIFT_MACRO_(_1, _2, _3, __VA_ARGS__, 0)
#define SHIFT_MACRO_(_1, _2, _3, NAME, ...) NAME
/* The extra macro adds an extra arguments to allow the smaller case. */

/**
 * Start a for-loop over the vertices of the graph or a subgraph.
 * @param idx Loop index; takes vertex values
 * @param subgraph (optional) The subgraph to loop over
 */
#define for_vertices(...) SHIFT_MACRO( \
		0, __VA_ARGS__, _for_vertices_of, _for_all_vertices)(__VA_ARGS__)

#define _for_all_vertices(idx) \
	for (vertex (idx) = 0; (idx) < N; ++(idx))

#define _for_vertices_of(idx, subgraph) \
	_for_all_vertices(idx) \
		if ((subgraph)[(idx)])

/**
 * Start a for-loop over pairs of vertices of the graph or a subgraph.
 * @param idx1
 * @param idx2
 * @param subgraph (optional)
 */
#define for_vertex_pairs(...) SHIFT_MACRO( \
		__VA_ARGS__, _for_vertex_pairs_of, _for_all_vertex_pairs)(__VA_ARGS__)

#define _for_all_vertex_pairs(idx1, idx2) \
	for (vertex (idx1) = 0; (idx1) < N; ++(idx1)) \
		for (vertex (idx2) = (idx1) + 1; (idx2) < N; ++(idx2))

#define _for_vertex_pairs_of(idx1, idx2, subgraph) \
	for (vertex (idx1) = 0; (idx1) < N; ++(idx1)) \
		if (contains(subgraph, idx1)) \
			for (vertex (idx2) = (idx1) + 1; (idx2) < N; ++(idx2)) \
				if (contains(subgraph, idx2))

/**
 * Start a for-loop over the edges the graph or a subgraph.
 * @param idx1
 * @param idx2
 * @param subgraph (optional)
 */
#define for_edges(...) SHIFT_MACRO( \
		__VA_ARGS__, _for_edges_of, _for_all_edges)(__VA_ARGS__)

#define _for_all_edges(idx1, idx2) \
	_for_all_vertex_pairs(idx1, idx2) \
		if (adj(idx1, idx2))

#define _for_edges_of(idx1, idx2, subgraph) \
	_for_vertex_pairs_of(idx1, idx2, subgraph) \
		if (adj(idx1, idx2))

/**
 * Start a for-loop over the neighbors of vertex u in the graph or a subgraph.
 * @param idx Index for neighbor loop
 * @param u Vertex whose neighbors are being used
 * @param subgraph (optional)
 * This will work regardless of whether the vertex is in the subgraph.
 */
#define for_neighbors(...) SHIFT_MACRO( \
		__VA_ARGS__, _for_neighbors_of, _for_all_neighbors)(__VA_ARGS__)

#define _for_all_neighbors(idx, u) \
	for (vertex (idx) = 0; (idx) < N; ++(idx)) \
		if (adj(u, idx))

#define _for_neighbors_of(idx, u, subgraph) \
	for (vertex (idx) = 0; (idx) < N; ++(idx)) \
		if (contains(subgraph, idx) && adj(u, idx))

/**
 * Check whether two vertices are adjacent in the graph.
 */
static inline bool adj(vertex u, vertex v)
{
	return ADJ[u][v];
}

/**
 * Check whether a subgraph contains a vertex
 */
static inline bool contains(const subgraph subgraph, const vertex u)
{
	return subgraph[u];
}

/**
 * Connect or disconnect two vertices of the graph
 */
static inline void set_edge(vertex u, vertex v, bool is_edge)
{
	ADJ[u][v] = ADJ[v][u] = is_edge;
}

#if N < 2
#error N too small, must be at least 2
#endif
#endif  /* GRAPH_H */
