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

## Parallelism

Since every iteration is an independent attempt - it only reads the
input graph, and doesn't depend on any other iteration's result -
`hpart` runs its iterations in parallel across threads using OpenMP,
rather than one at a time on a single core. Each thread repeatedly
builds and locally improves its own partition for the full run-time
budget; the only thing threads share is the read-only input graph, so
they never interfere with each other. Once every thread's time budget
is up, the lowest-cost partition found by *any* thread is reported.

This doesn't change what each individual iteration does - it just buys
more of them in the same wall-clock time. With N threads, roughly N
times as many independent attempts fit in the same run-time budget,
which (per the thesis's own finding above) tends to produce a better
result, since it's equivalent to giving the algorithm proportionally
more running time. Results are consistent but not bit-for-bit
reproducible across different thread counts: more threads mean more
total attempts, and the specific partitions explored depend on the
RNG stream each thread happens to run, not just the single global seed
used.

The number of threads defaults to the number of available cores, or
can be capped with the standard `OMP_NUM_THREADS` environment
variable, e.g. `OMP_NUM_THREADS=4 ./hpart ...`.

## Files

- `hmake.c` / `hmake.h` — generates a random 0-1 graph and writes it in
  the input format `hpart` expects.
- `hpart.c` — CLI argument parsing and the driver loop, parallelized
  across threads with OpenMP (see "Parallelism" above).
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

`hpart` requires a compiler and standard library with OpenMP support
(e.g. GCC or Clang with `libgomp`/`libomp` installed) - the Makefile
passes `-fopenmp` when building it.

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
- `run-time` — wall-clock seconds to search before reporting the best
  partition found, across all threads (see "Parallelism" above).
