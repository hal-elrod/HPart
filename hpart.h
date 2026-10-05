/* .........................................................................
   hpart.h: shared types, constants, globals, and function prototypes
   for hpart.c (CLI + driver loop), readpart.c (graph input), and
   greedy.c (partition construction + local-search swap strategies).

   This used to be one file, HPART.C. The original (1989) toolchain for
   this program ran on hardware that didn't support splitting code
   across .h/.c files at all, so everything lived in a single file.
   ....................................................................... */

#ifndef HPART_H
#define HPART_H

#include <stdio.h>

#define ind(i,j,nn) (i*nn+j)

/* Per-thread RNG stream for the OMP-parallel attempt loop in main().
   Each thread sets its own rng_seed once (derived from the single
   RNG_SEED below plus its thread number) and then only ever touches
   its own copy via rand_r(), so construction randomness stays both
   race-free and reproducible instead of every thread racing on libc's
   global rand()/srand() state. */
extern __thread unsigned int rng_seed;
#undef random
#define random(num) ((int)(rand_r(&rng_seed) % (num)))

/* Upper bound on cand_list_size (the <candidate-list-size> command-line
   argument). greedypart's cand[]/candcost[] and heappart's clist[]/
   templist[] are all sized off this constant - keep it in sync with
   those array sizes if they ever change, or candidate-list writes can
   run off the end. */
#define MAX_CAND_LIST_SIZE 9

/* "Infinity" sentinels. A node's gain (costa[]/heapslot.alpha) is
   bounded by its degree, which can be at most nn-1 - these constants
   need to stay far above that for every nn this program can practically
   run on, not just the small graphs it was originally tested with. They
   are only ever compared against, never added together, so there's no
   overflow risk in giving them this much headroom within int's range.
   Used to seed a min-search (with POS_INF - anything real is lower) or
   a max-search (with NEG_INF - anything real is higher), and to mark a
   candidate slot as "not holding a real value yet". */
#define POS_INF 2000000000
#define NEG_INF -2000000000

/* Margin below/above the true +-POS_INF/NEG_INF sentinel used to tell
   whether a heap slot's gain (heapslot.alpha) is still a real, active
   value or has been marked removed (set to exactly -POS_INF/POS_INF).
   Kept as a separate constant with headroom rather than comparing
   against the sentinel exactly: a removed slot's alpha keeps getting
   nudged by +-1 as its neighbors are placed (see heappart's update
   loops), so the margin needs to comfortably exceed the largest degree
   this program can practically handle. */
#define REMOVED_SLOT_MARGIN 1000000000

/* Initial "worst case" gain in aslightswap/slightswap/slightestswap's
   search for the smallest positive-gain swap: larger than any gain a
   real swap could produce (bounded by node degree), so the first
   positive-gain swap found always replaces it. */
#define WORST_SWAP_SENTINEL 1000000000

/* Initial value of main()'s mincval, the best (lowest) partition cost
   seen so far: larger than any real partition cost (bounded by the
   number of edges), so the first attempt's result always replaces it.
   If this is ever printed as the final "min cost", no attempt beat it -
   a sign this needs to be raised further for the graph size in use. */
#define INITIAL_MIN_COST 2000000000

/* Fixed RNG seed so repeated runs over the same inputs and parameters
   are reproducible, which matters for comparing partitioning
   strategies against each other. */
#define RNG_SEED 32063

/* Defined in hpart.c (set there from parse_args()'s result), shared
   with readpart.c (getgraph/readgraph) and greedy.c (greedypart/
   heappart). cand_list_size is the validated <candidate-list-size>
   argument; big_flag is 1 unless modeb is 4 (the compact, no-matrix
   swap strategy), in which case the igraph[] adjacency matrix is
   never allocated or populated. */
extern int cand_list_size;
extern int big_flag;

/* Parsed, validated command-line arguments. parse_args() exits the
   program (after printing usage/an error) rather than returning an
   invalid value, so every field here is guaranteed valid on return. */
typedef struct {
	char *inputfile;
	int modea;		/* 1 or 2 */
	int modeb;		/* 1-4 */
	int cand_list_size;	/* 1..MAX_CAND_LIST_SIZE */
	float run_time;		/* > 0 */
	} cmdargs;

/* Parses and validates argv, printing usage and exiting on any problem
   (including argc being wrong) instead of returning. */
cmdargs parse_args(int argc, char *argv[]);

typedef struct anode{
	int node;
	struct anode *next;
	}nodez;

typedef nodez *nodep;

/* greedypart's sindex[]: a doubly-linked free list, indexed by node id,
   of the nodes not yet placed into set A or B. */
typedef struct {
	int back;
	int next;
	} linknode;

/* heappart's heapa[]/heapb[]: one slot of a binary heap. hnode is the
   graph node occupying this slot; alpha is its current gain. heapa is
   max-on-top (best node to add to set A), heapb is min-on-top. */
typedef struct {
	int hnode;
	int alpha;
	} heapslot;

/* heappart's index[]: for a given graph node, where it currently sits
   in each heap, so a node's gain can be updated in both heaps in O(log n)
   without searching. apoint = slot in heapa, bpoint = slot in heapb. */
typedef struct {
	int apoint;
	int bpoint;
	} crossindex;

/* Opens inputfile and allocates the shared, read-only graph data:
   igraph[] (if big_flag) and alist[]. nn/ne are filled in from the
   file's header line. */
void getgraph(char *inputfile,int **igraph,int *nn,int *ne,nodez **alist);

/* Allocates one thread's private working arrays for building and
   holding a partition: ma[]/mb[] (the two sides) and sindex[]
   (greedypart's free-list of unplaced nodes). Call once per thread and
   reuse across that thread's attempts, same as the old single-threaded
   main() did once for the whole program. */
void alloc_partition(int nn,int **ma,int **mb,linknode **sindex);

/* Fills igraph (if big_flag) and alist from the input file's edge list. */
void readgraph(int ne,int nn,int igraph[],nodez alist[]);

/* Builds an initial 2-partition greedily, picking each node from a
   randomized candidate list, using a doubly-linked list of unplaced
   nodes (sindex) to find candidates in O(n) per pass. */
void greedypart(int costa[],int nn,int ma[],
		int mb[],linknode sindex[],nodez alist[]);

/* Builds an initial 2-partition greedily, same as greedypart but using
   two gain heaps (heapa/heapb) instead of a linear scan to pick each
   next candidate. */
void heappart(int costa[],int nn,int ma[],int mb[],nodez alist[]);

/* Returns elapsed wall-clock seconds since the previous call. */
float periodt(void);

/* Allocates and zeroes the per-attempt gain array costa[]. */
void remem(int **pa,int nn);

/* Local-search postprocessors: repeatedly swap one node from set A with
   one from set B whenever doing so improves the cross-edge weight,
   until no improving swap remains. They differ only in which improving
   swap they pick each pass (and in how edge weight is looked up):
     aslightswap    - smallest gain under 3, via the adjacency list
                       (used when big_flag is off / mode 4)
     hswap          - first improving swap found, via igraph[]
     slightswap     - smallest gain under 3, via igraph[]
     slightestswap  - smallest gain exactly 1, via igraph[] */
void aslightswap(int ma[],int mb[],int costa[],int nn,int *cval,nodez alist[]);
void slightswap(int igraph[],int ma[],int mb[],int costa[],
		int nn,int *cval,nodez alist[]);
void slightestswap(int igraph[],int ma[],int mb[],int costa[],
		   int nn,int *cval,nodez alist[]);
void hswap(int igraph[],int ma[],int mb[],int costa[],int nn,
	   int *cval,nodez alist[]);

#endif
