/* .........................................................................
   HPART.h: shared types and function prototypes for HPART.c

   Split out of HPART.c for clarity. The original (1989) toolchain for this
   program ran on hardware that didn't support splitting code across .h/.c
   files, so everything used to live in one file.
   ....................................................................... */

#ifndef HPART_H
#define HPART_H

#include <stdio.h>

#define ind(i,j,nn) (i*nn+j)
#undef random
#define random(num) (rand() % (num))

/* Upper bound on cand_list_size (the <cl-size> command-line argument).
   greedypart's cand[]/candcost[] and heappart's clist[]/templist[] are
   all sized off this constant - keep it in sync with those array sizes
   if they ever change, or candidate-list writes can run off the end. */
#define MAX_CAND_LIST_SIZE 9

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

/* Opens inputfile and allocates the graph/working arrays. */
void getgraph(char *inputfile,
	      int  **igraph,int **a, int **b,linknode **sindex,
	      int *nn,int *ne,nodez **alist);

/* Fills igraph (if big_flag) and alist from the input file's edge list. */
void readgraph(int ne,int nn,int igraph[],nodez alist[]);

/* Builds an initial 2-partition greedily, picking each node from a
   randomized candidate list, using a doubly-linked list of unplaced
   nodes (sindex) to find candidates in O(n) per pass. */
void greedypart(int costa[],int nn,int ma[],
		int mb[],linknode sindex[],nodez alist[]);

/* Sift a slot toward the root of heapa (max-on-top) / heapb (min-on-top)
   after its gain has increased / decreased. */
void upaheap(heapslot heap[],int i,crossindex index[]);
void upbheap(heapslot heap[],int i,crossindex index[]);

/* Sift a slot toward the leaves of heapa (max-on-top) / heapb (min-on-top)
   after its gain has decreased / increased. */
void downaheap(heapslot heap[],int i,int nn,crossindex index[]);
void downbheap(heapslot heap[],int i,int nn,crossindex index[]);

/* Builds an initial 2-partition greedily, same as greedypart but using
   two gain heaps (heapa/heapb) instead of a linear scan to pick each
   next candidate. */
void heappart(int costa[],int nn,int ma[],int mb[],nodez alist[]);

/* Returns elapsed wall-clock seconds since the previous call. */
float periodt(void);

/* Allocates and zeroes the per-attempt gain array costa[]. */
void remem(int **pa,int nn);

/* Returns 1 if node1 and node2 are adjacent in the input graph, else 0.
   Used instead of the igraph[] adjacency matrix when big_flag is off. */
int finda_weight(nodez alist[],int node1,int node2);

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
