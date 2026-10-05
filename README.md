# A GRASP for the 2-Partition Problem

An implementation of a Greedy Randomized Adaptive Search Procedure
(GRASP) for the network 2-partition problem, from Hal Elrod's 1989
M.S. thesis at the University of Texas at Austin, *A Greedy Randomized
Adaptive Search Procedure for the 2-partition Problem* (included here
as `H. Elrod THESIS.pdf`).

## The problem

Given a graph with weighted edges, divide its nodes into two
equal-sized sets so that the total weight of the edges crossing
between the two sets is as small as possible. This is NP-hard, and
comes up whenever something needs to be physically or logically split
into two halves while minimizing the connections between them — e.g.
placing VLSI circuit modules on a chip, or partitioning a program's
subroutines across memory pages to reduce paging between them.

## The approach

Each GRASP iteration has two phases:

1. **Construction** — build a partition one node at a time. At each
   step, instead of always picking the single best-looking node
   (purely greedy, deterministic), a short candidate list of the best
   few nodes is formed and one is picked at random. This is what makes
   different iterations explore different parts of the solution space.
2. **Local search** — repeatedly swap a node from one side with a node
   from the other whenever doing so lowers the cut weight, until no
   such swap is left.

Many iterations are run within a time budget, and the best partition
found across all of them is reported. The thesis compares this
approach against the then-dominant Kernighan-Lin (K&L) heuristic
across thousands of random and geometric graph instances, and finds
GRASP wins on the large majority of them — more decisively the more
running time it's given, since each additional iteration is an
independent shot at a better local optimum.

## Files

- `hmake.c` / `hmake.h` — generates a random 0-1 graph and writes it in
  the input format `hpart` expects.
- `hpart.c` — CLI argument parsing and the driver loop.
- `readpart.c` — reads the graph file (`getgraph`/`readgraph`).
- `greedy.c` — partition construction (`greedypart`/`heappart`) and the
  local-search swap strategies (`aslightswap`/`hswap`/`slightswap`/
  `slightestswap`).
- `hpart.h` — shared types, constants, and prototypes for the three
  `hpart` source files.
- `Makefile` — builds both `hmake` and `hpart`.
- `H. Elrod THESIS.pdf` — the thesis this code is from; see it for the
  full algorithm description, complexity analysis, and experimental
  results.

## Building

```
make
```

builds `hmake` and `hpart`. `make clean` removes the binaries and
object files.

## Usage

First generate a graph:

```
./hmake graph.txt
```

It will prompt for a random seed, a node count, and an edge
probability, then write the graph to `graph.txt`.

Then partition it:

```
./hpart graph.txt <modea> <modeb> <candidate-list-size> <run-time>
```

Run `./hpart` with no arguments (or `-h`/`--help`) for a full
explanation of each argument. Briefly:

- `modea` — initial construction method: `1` heap-greedy, `2` greedy
  (linear scan).
- `modeb` — local-search swap strategy: `1` first-improvement
  2-exchange, `2` slight swap, `3` slightest swap, `4` compact slight
  swap (adjacency-list based, avoids the O(n²) adjacency matrix — use
  for very large/sparse graphs).
- `candidate-list-size` — how many top candidates the construction
  phase considers at each step before picking one at random (1-9). The
  thesis found 2-4 gives the best balance between randomness and
  solution quality; much higher values degrade the initial solutions
  passed to local search.
- `run-time` — seconds to search before reporting the best partition
  found.
