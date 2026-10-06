 /* .........................................................................
	HPART.C: A GRASP approach to the 2-partition problem

	GRASP == Greedy Randomized Adaptive Search Procedure
	Builds a low weight partition of 0-1 graph by greedily adding
	pairs of nodes from a candidate list of those nodes the maximize
	the current partition, then the weight of the partition is reduced
	by exchanging pairs of nodes when the exchange will increase the
	weight of the inner edges.

	This file holds CLI argument parsing and the driver loop. Graph
	input lives in readpart.c; partition construction and the
	local-search swap strategies live in greedy.c.

	Hal Elrod   --  Operations Research Group
			Department of Mechanical Engineering
			The University of Texas at Austin
			Austin, TX  78712

	Spring Break, March 1989
	last modified 6/30/89
   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   Input:	number of nodes, number of edges
		node, node, weight
		 .      .     .
		 .	.     .
		node, node, weight

   Use:		HPART <inputfile> <modea> <modeb> <c-list> <run-time>
   Where:	<modea> = "1", heap greedy unmatched partition
			= "2", greedy unmatched partition
		<modeb> = "1", first swap -- generic 2-exchange
			= "2", slight swap
			= "3", slightest swap
			= "4", compact slight swap

   ....................................................................... */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <pthread.h>
#include "hpart.h"

static clock_t oldmtime;
int cand_list_size,big_flag = 1;
__thread unsigned int rng_seed;
/* Guards combining each thread's local_mincval/local_attemps into the
   shared mincval/num_attemp at the end of main()'s parallel region. A
   plain pthread_mutex instead of #pragma omp critical, because GCC's
   libgomp implements omp critical/reduction in a way ThreadSanitizer
   doesn't fully recognize as synchronization (confirmed false-positive,
   not a real race, but a real mutex removes the ambiguity entirely). */
static pthread_mutex_t mincval_mutex = PTHREAD_MUTEX_INITIALIZER;

/* Prints a full explanation of every argument. Shown on a missing/bad
   argument, or on request via -h/--help/-?. */
static void print_usage(const char *prog)
{
	printf("Usage: %s <inputfile> <modea> <modeb> <candidate-list-size> <run-time>\n\n",prog);
	puts("  inputfile             Graph file as written by hmake: first line");
	puts("                        \"n,ne\", then ne lines of \"node,node,weight\".");
	puts("");
	puts("  modea                 Initial partition construction method:");
	puts("                          1  heap-greedy (gain heaps; faster candidate lookup)");
	puts("                          2  greedy (linear scan for each candidate)");
	puts("");
	puts("  modeb                 Local-search swap strategy used to improve the");
	puts("                        partition:");
	puts("                          1  first-improvement 2-exchange");
	puts("                          2  slight swap (best swap with gain < 3)");
	puts("                          3  slightest swap (best swap with gain == 1)");
	puts("                          4  compact slight swap (adjacency-list based,");
	puts("                             no n x n matrix - use for very large/sparse");
	puts("                             graphs)");
	puts("");
	printf("  candidate-list-size   Candidate list size, 1-%d. Larger values consider\n",
	       MAX_CAND_LIST_SIZE);
	puts("                        more candidates per step at some cost in speed.");
	puts("");
	puts("  run-time              Seconds to search before reporting the best");
	puts("                        partition found (the program keeps looking for");
	puts("                        roughly 2x this long before it actually stops).");
	puts("");
	printf("Example:\n  %s graph.txt 1 1 %d 5\n",prog,MAX_CAND_LIST_SIZE);
}

/* Parses and validates argv. On any problem - wrong argument count, or
   any argument out of range - prints an explanation and exits rather
   than returning, so every field of the result is guaranteed usable. */
cmdargs parse_args(int argc, char *argv[])
{
cmdargs args;
char *end;
long cl;

	if (argc == 2 && (strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0
			  || strcmp(argv[1],"-?") == 0))
		{
		print_usage(argv[0]);
		exit(0);
		}

	if (argc != 6)
		{
		if (argc > 1)
			printf("Error: expected 5 arguments, got %d.\n\n",argc - 1);
		print_usage(argv[0]);
		exit(argc > 1 ? 1 : 0);
		}

	args.inputfile = argv[1];

	if (strlen(argv[2]) != 1 || (argv[2][0] != '1' && argv[2][0] != '2'))
		{
		printf("Error: modea must be 1 or 2, got \"%s\".\n\n",argv[2]);
		print_usage(argv[0]);
		exit(1);
		}
	args.modea = argv[2][0] - '0';

	if (strlen(argv[3]) != 1 || argv[3][0] < '1' || argv[3][0] > '4')
		{
		printf("Error: modeb must be 1, 2, 3, or 4, got \"%s\".\n\n",argv[3]);
		print_usage(argv[0]);
		exit(1);
		}
	args.modeb = argv[3][0] - '0';

	cl = strtol(argv[4],&end,10);
	if (*end != '\0' || end == argv[4] || cl < 1 || cl > MAX_CAND_LIST_SIZE)
		{
		printf("Error: candidate-list-size must be an integer from 1 to %d, got \"%s\".\n\n",
		       MAX_CAND_LIST_SIZE,argv[4]);
		print_usage(argv[0]);
		exit(1);
		}
	args.cand_list_size = (int)cl;

	args.run_time = strtof(argv[5],&end);
	if (*end != '\0' || end == argv[5] || args.run_time <= 0)
		{
		printf("Error: run-time must be a positive number of seconds, got \"%s\".\n\n",
		       argv[5]);
		print_usage(argv[0]);
		exit(1);
		}

	return args;
}

/* Returns elapsed wall-clock seconds since the previous call to periodt(). */
float periodt()
{
 clock_t marktime,temp;

	temp = oldmtime;
	marktime = clock();
	oldmtime = marktime;
	return	((marktime-temp)/CLOCKS_PER_SEC);
}

/* Allocates and zeroes the per-attempt gain array costa[]. */
void remem(int **pa,int nn)
{
	*pa = (int *) calloc (nn+1,sizeof(int));
	if (*pa == NULL)
		{
		puts("Couldn't allocate an array");
		exit(0);
		}
}

/* Driver: reads the graph once, then runs one attempt loop per thread
   in parallel - each building a fresh partition (modea) and locally
   improving it (modeb), for run_time*2 seconds of wall-clock time -
   and reports the lowest-cost partition seen across all of them.

   Threads share the read-only graph (igraph[]/alist[]) but each keeps
   its own partition (ma[]/mb[]/sindex[]), gain array (costa[]), and RNG
   stream (rng_seed), so attempts never interfere with each other. */
int main (int argc,char *argv[])
{
int nn,ne,num_attemp = 0,mincval = INITIAL_MIN_COST;
int *igraph;
nodez *alist;
cmdargs args;
double run_time,start_time;

	args = parse_args(argc,argv);
	cand_list_size = args.cand_list_size;
	if(args.modeb == 4)
		big_flag = 0;
	run_time = args.run_time;

	getgraph(args.inputfile,&igraph,&nn,&ne,&alist);
	readgraph(ne,nn,igraph,alist);

	start_time = omp_get_wtime();

	#pragma omp parallel default(none) \
		shared(args,igraph,alist,nn,run_time,start_time,mincval,num_attemp,mincval_mutex)
		{
		int *ma,*mb,*costa,cval;
		linknode *sindex;
		int local_mincval = INITIAL_MIN_COST,local_attemps = 0;

		rng_seed = RNG_SEED + (unsigned int)omp_get_thread_num();
		alloc_partition(nn,&ma,&mb,&sindex);

		while (omp_get_wtime() - start_time < run_time * 2)
			{
			local_attemps++;
			remem(&costa,nn);
			switch (args.modea)
				{
				case 1: heappart(costa,nn,ma,mb,alist);
					break;
				case 2: greedypart(costa,nn,ma,mb,sindex,alist);
					break;
				}
			cval = 0;
			switch (args.modeb)
				{
				case 1: hswap(igraph,ma,mb,costa,nn,&cval,alist);
					break;
				case 2: slightswap(igraph,ma,mb,costa,nn,&cval,alist);
					break;
				case 3: slightestswap(igraph,ma,mb,costa,nn,&cval,alist);
					break;
				case 4: aslightswap(ma,mb,costa,nn,&cval,alist);
					break;
				}
			free(costa);
			if (cval < local_mincval)
				local_mincval = cval;
			}    /* while */
		free (ma); free(mb);
		free (sindex);

		/* Combine this thread's results into the shared totals. */
		pthread_mutex_lock(&mincval_mutex);
		if (local_mincval < mincval)
			mincval = local_mincval;
		num_attemp += local_attemps;
		pthread_mutex_unlock(&mincval_mutex);
		}    /* omp parallel */

	printf("min cost = %d\n",mincval);
	printf("number of attempts = %d\n",num_attemp);
	freegraph(nn,alist);
	free (alist);
	if (big_flag) free(igraph);
	return EXIT_SUCCESS;
}
