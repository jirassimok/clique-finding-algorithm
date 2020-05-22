/**
 * @file clique.c
 * @author Jacob Komissar
 * @date 2020-05-02
 *
 * Find cliques in a graph with a planted clique.
 *
 * Compile this program with 'gcc -Wall -Wextra -O3 --std=c99 clique.c'
 *
 * In this program, the graph has constant size and is represented as an
 * adjacency matrix.
 *
 * Subgraphs are represented as boolean arrays indicating whether each vertex
 * is in the subgraph.
 */
// This program uses static variables---don't even think about multi-threading
#include <limits.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "graph.h"

#define CAT_(x, y, z) x##y##_##z
#define CAT(x, y, z) CAT_(x, y, z)
#define MACRO_ID(x) CAT(_macro_generated_id_, x, __LINE__)

#define dotimes(n) for (int MACRO_ID(n) = 0; MACRO_ID(n) < n; ++MACRO_ID(n))


/**
 * Default maximum size of planted clique
 */
static const int MAX_K = 128;

/**
 * Default minimum size of planted clique
 */
static const int MIN_K = 28;

/**
 * Default number of trials to run for each clique size
 */
static const int TRIALS_PER_K = 1000;


void randomize_graph(double);
void plant_clique(int);
int rand_under(int);

int ldr(void);

bool is_clique(const subgraph);
int size(const subgraph);
int min_deg(const subgraph, int*);
void reset_degrees(const subgraph, int*);

void add_vertex(vertex, subgraph, int*);
void remove_vertex(vertex, subgraph, int*);

int deg(vertex, const subgraph);

#define READ_INT_ARG(idx, name) \
	if (argc > idx) \
		if (sscanf(argv[idx], "%d", &name) < 1) \
			fprintf(stderr, "Using default value %d for %s\n", name, #name);

int main(int argc, char *argv[])
{
	int trials = TRIALS_PER_K;
	int min_k = MIN_K;
	int max_k = MAX_K;
	READ_INT_ARG(1, trials);
	READ_INT_ARG(2, min_k);
	READ_INT_ARG(3, max_k);

	fprintf(stderr, "Config: %d trials for each k from %d to %d\n", trials, min_k, max_k);

	struct timeval now;
	gettimeofday(&now, NULL);
	srand48(now.tv_sec * 1000000 + now.tv_usec);

	double means[max_k + 1];

	for (int k = min_k; k <= max_k; ++k) {
		fprintf(stderr, "Running with %d-cliques.", k);

		int total_size = 0;
		for (int i = 1; i <= trials; ++i) {
			if (i % 10 == 0) {
				fprintf(stderr, "\rRunning with %d-cliques. %d ", k, i);
			}

			randomize_graph(0.5);
			plant_clique(k);
			total_size += ldr();
		}

		means[k] = (double)total_size / (double)trials;
		fprintf(stderr, "\rRunning with %d-cliques. Mean found: %f\n", k, means[k]);
	}

	// Print to stdout to save in a convenient format
	printf("# Config: %d trials for each k from %d to %d\n", trials, min_k, max_k);
	for (int k = min_k; k <= max_k; ++k) {
		printf("%d %f\n", k, means[k]);
	}

	return 0;
}

/**
 * Randomly set each edge of the graph with probability p.
 */
void randomize_graph(double p)
{
	for_vertex_pairs(u, v) {
		set_edge(u, v, drand48() < p);
	}
}

/**
 * Plant a random size-k clique in the global graph.
 */
void plant_clique(int k)
{
	static int indices[N]; // vertices indices, in no particular order
	int vertices[k]; // vertices in the clique; a proper index array
	// TODO: the vertices array isn't necessary; after the shuffle it's the first k elements of the other array

	if (!indices[1]) {
		// Initialize indices to {0, 1, ..., N-1}
		// This only ever must be done once; the list can be in any order
		for (int i = 0; i < N; ++i) {
			indices[i] = i;
		}
	}

	// Choose vertices for the clique using partial Fisher-Yates shuffle
	for (int i = 0; i < k; ++i) {
		int j = rand_under(N - i);
		vertex u = indices[i + j];

		vertices[i] = u; // save for later use

		indices[i + j] = indices[i];
		indices[i] = u;
	}

	// Make the chosen vertices into a clique
	for (int i = 0; i < k; ++i) {
		vertex u = vertices[i];
		for (int j = i + 1; j < k; ++j) {
			vertex v = vertices[j];
			set_edge(u, v, true);
		}
	}
}

/**
 * Generate a random integer on [0, limit).
 */
int rand_under(int limit)
{
	const long rand_limit = 1L << 31;  // lrand48 generates on [0, 2^31 - 1]

	long bias_splits = rand_limit / limit;
	long bias_limit = bias_splits * limit;

	// If the random number doesn't fit in an even multiple of the limit, try
	// again. This is largely efficient, with worst case expected tries
	// slightly below 2 when the limit is just under half of rand_limit.
	int rand;
	do {
		rand = (int)lrand48();
	} while (rand >= bias_limit);

	return rand / bias_splits;
}

/**
 * Run the low-degree-removal algorithm on the global graph to find a clique.
 *
 * @return The size of the located clique.
 */
int ldr()
{
	static subgraph V;
	static subgraph G;
	// Make V empty and G full
	memset(V, false, N * sizeof (*V));
	memset(G, true, N * sizeof (*G));

	// Store the degrees in this array to save time.
	static int deg[N];
	reset_degrees(G, deg);

	// Removal phase
	while (!is_clique(G)) {
		vertex u = min_deg(G, deg);
		add_vertex(u, V, NULL);
		remove_vertex(u, G, deg);
	}

	// Inclusion phase
	reset_degrees(G, deg);
	int size_G = size(G);
	for_vertices(v, V) {
		// v is adjavent to all of G iff its degree in G is as large as G
		if (deg[v] == size_G) {
			size_G += 1;
			add_vertex(v, G, deg);
		}
	}

	return size_G;
}

/**
 * Write to the given array the degrees in the given graph of all vertices.
 */
void reset_degrees(const subgraph G, int *deg)
{
	memset(deg, 0, N * sizeof (int));
	for_vertices(u, G) {
		for_neighbors(v, u) {
			deg[v] += 1;
		}
	}
}

/**
 * Add the given vertex to the given subgraph, updating the given degrees if
 * not NULL.
 */
void add_vertex(vertex u, subgraph G, int *deg)
{
	G[u] = true;
	if (deg != NULL) {
		for_neighbors(v, u) {
			deg[v] += 1;
		}
	}
}

/**
 * Remove the given vertex from the given subgraph, updating the given degrees
 * if not NULL.
 */
void remove_vertex(vertex u, subgraph G, int *deg)
{
	G[u] = false;
	if (deg != NULL) {
		for_neighbors(v, u) {
			deg[v] -= 1;
		}
	}
}

/**
 * Find the index of the vertex of minimum degree in the given subgraph.
 */
int min_deg(const subgraph G, int *deg)
{
	int d_min = -1;
	vertex u_min = -1;
	for_vertices(u, G) {
		int d = deg[u];
		if (d < d_min || d_min < 0) {
			u_min = u;
			d_min = d;
		}
	}
	return u_min;
}

/**
 * Get the degree of the given vertex in the given subgraph.
 */
int deg(vertex u, const subgraph G)
{
	int d = 0;
	for_neighbors(v, u, G) {
		d += 1;
	}
	return d;
}


/**
 * Check if a subgraph of the global graph is a clique
 */
bool is_clique(const subgraph G)
{
	for_vertex_pairs(u, v, G) {
		if (!adj(u, v)) {
			return false;
		}
	}
	return true;
}

/**
 * Get the number of vertices in a subgraph.
 */
int size(const subgraph G)
{
	int n = 0;
	for_vertices(u, G) {
		n += 1;
	}
	return n;
}
