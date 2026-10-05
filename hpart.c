 /* .........................................................................
	HPART.C: A GRASP approach to the 2-partition problem

	GRASP == Greedy Randomized Adaptive Search Procedure
	Builds a low weight partition of 0-1 graph by greedily adding
	pairs of nodes from a candidate list of those nodes the maximize
	the current partition, then the weight of the partition is reduced
	by exchanging pairs of nodes when the exchange will increase the
	weight of the inner edges.

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

   For MS-DOS systems, HPART should be compiled under COMPACT or HUGE models.
   ....................................................................... */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hpart.h"

/* CLK_TCK was removed from modern glibc headers; it was always equal
   to CLOCKS_PER_SEC on the systems this ran on. */
#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

static FILE *f_in;
clock_t oldmtime;
int cand_list_size,big_flag = 1;

/* Prints a full explanation of every argument. Shown on a missing/bad
   argument, or on request via -h/--help/-?. */
static void print_usage(const char *prog)
{
	printf("Usage: %s <inputfile> <modea> <modeb> <cl-size> <run-time>\n\n",prog);
	puts("  inputfile  Graph file as written by hmake: first line \"n,ne\",");
	puts("             then ne lines of \"node,node,weight\".");
	puts("");
	puts("  modea      Initial partition construction method:");
	puts("               1  heap-greedy (gain heaps; faster candidate lookup)");
	puts("               2  greedy (linear scan for each candidate)");
	puts("");
	puts("  modeb      Local-search swap strategy used to improve the partition:");
	puts("               1  first-improvement 2-exchange");
	puts("               2  slight swap (best swap with gain < 3)");
	puts("               3  slightest swap (best swap with gain == 1)");
	puts("               4  compact slight swap (adjacency-list based, no");
	puts("                  n x n matrix - use for very large/sparse graphs)");
	puts("");
	printf("  cl-size    Candidate list size, 1-%d. Larger values consider more\n",
	       MAX_CAND_LIST_SIZE);
	puts("             candidates per step at some cost in speed.");
	puts("");
	puts("  run-time   Seconds to search before reporting the best partition");
	puts("             found (the program keeps looking for roughly 2x this");
	puts("             long before it actually stops).");
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
		printf("Error: cl-size must be an integer from 1 to %d, got \"%s\".\n\n",
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

/* Driver: reads the graph once, then repeatedly builds a fresh partition
   (modea) and locally improves it (modeb) for run_time*2 seconds,
   keeping track of the lowest-cost partition seen, then reports it.
   (The loop body's "if (!oneflag && ...)" block for printing interim
   stats at the run_time mark never runs - oneflag is initialized to 1,
   so !oneflag is always false - this appears to be leftover from an
   earlier version of the instrumentation.) */
int main (int argc,char *argv[])
{
double sigx = 0.0 ,sigx2 = 0.0, sig;
int nn,ne,cval = 0,num_attemp = 0,mincval = 32600;
int *igraph, *ma,*mb,*costa,oneflag = 1;  /* ,*check,x; */
nodez *alist;
linknode *sindex;
clock_t startt;
cmdargs args;
float matcht = 0,partt= 0,swapt = 0,bestt,ttime,curtime,run_time;

	args = parse_args(argc,argv);
	cand_list_size = args.cand_list_size;
	if(args.modeb == 4)
		big_flag = 0;
	run_time = args.run_time;

	periodt();
	startt = clock();
	srand(32063);
	getgraph(args.inputfile,&igraph,&ma,&mb,&sindex,
		 &nn,&ne,&alist);

	readgraph(ne,nn,igraph,alist);
	curtime = clock()/CLK_TCK;
	while (curtime < run_time * 2)
		{
		num_attemp++;
		remem(&costa,nn);
		switch (args.modea)
			{
			case 1: heappart(costa,nn,ma,mb,alist);
				partt += periodt();
				break;
			case 2: greedypart(costa,nn,ma,mb,sindex,alist);
				partt += periodt();
				break;
			}
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
		swapt += periodt();
		free(costa);
		sig = cval;
		sigx += sig;
		sigx2 += (sig * sig);
		if (cval < mincval)
			{
			bestt = ((clock() - startt)/CLK_TCK);
			mincval = cval;
			}
		cval = 0;
		curtime = clock()/CLK_TCK;
		if (!oneflag && curtime >= run_time)
			{
			oneflag++;
			ttime = ((clock() - startt)/CLK_TCK);
			printf("%d\n",mincval);
			printf("%f\n",ttime);
			printf("%f\n",bestt);
			printf("%f\n",matcht);
			printf("%f\n",partt);
			printf("%f\n",swapt);
			printf("%lg\n",sigx);
			printf("%lg\n",sigx2);
			printf("%d\n",num_attemp);
			}
		}    /* while */
        ttime = ((clock() - startt)/CLK_TCK);
	printf("min cost = %d\n",mincval);
	printf("number of attempts = %d\n",num_attemp);
	free (alist);
	if (big_flag) free(igraph);
	free (ma); free(mb);
	free (sindex);
	return EXIT_SUCCESS;
}

/*  Opens inputfile, inputs nn, ne, and allocates memory for matching. */
void getgraph(char *inputfile,
	      int **igraph,int **ma, int **mb, linknode **sindex,
	      int *nn,int *ne,nodez **alist)
{
int igraph_size,nn1,nnhalf;

	f_in = fopen(inputfile,"r");
	if (f_in == NULL)
		{
		printf("Error: couldn't open input file \"%s\".\n",inputfile);
		exit(1);
		}
	fscanf(f_in,"%d,%d",nn,ne);
        nn1 = (*nn)+1;
	nnhalf = (*nn) / 2;
	igraph_size = (nn1) * (nn1);
	if(big_flag)
		{
		*igraph = (int *) calloc (igraph_size,sizeof(int));
		if (*igraph == NULL)
			{
			puts("Couldn't make nxn array");
			exit(0);
			}
		}
	*ma = (int *) calloc (nnhalf + 1,sizeof(int));
	*mb = (int *) calloc (nnhalf + 1,sizeof(int));
	*sindex = (linknode *) calloc (nn1,sizeof(linknode));
	*alist = (nodez *) calloc (nn1,sizeof(nodez));
	if (*ma == NULL || *mb == NULL || *sindex == NULL || *alist == NULL)
		{
		puts("Couldn't allocate an array");
		exit(0);
		}
}

/* Reads the ne edges from the input file, building the adjacency list
   (alist) for every node and, if big_flag is set, the igraph[] adjacency
   matrix used for O(1) edge-weight lookups. */
void readgraph(int ne,int nn,int igraph[],nodez alist[])
{
int x,i,j,weight;
nodez *newi,*newj;
nodep *head;

	head = (nodep *) calloc (nn+1,sizeof(nodep));
	if(head == NULL)
		{
		puts("Couldn't allocate array in readgraph");
		exit(0);
		}
	for (x = 1; x<= nn; x++)
		{
		alist[x].next = NULL;
		head[x] = &alist[x];
		alist[x].node = 0;
		}

	for (x=1;x < ne + 1;x++)
		{
		fscanf(f_in,"%d,%d,%d",&i,&j,&weight);
		/* igraph[i]++; */
		/* igraph[j]++; */
		if (big_flag)
			{
			igraph[ind(i,j,nn)] = weight;
			igraph[ind(j,i,nn)] = weight;
			}
		alist[i].node++;
		alist[j].node++;
		newi = (nodez *) malloc (sizeof(nodez));
		newj = (nodez *) malloc (sizeof(nodez));
		if (newi == NULL || newj == NULL)
			{
			puts("Couldn't allocate node in readgraph");
			exit(0);
			}
		head[i]->next = newi;
		head[j]->next = newj;
		newi->node = j;
		newj->node = i;
		newi->next = NULL;
		newj->next = NULL;
		head[i] = newi;
		head[j] = newj;
		}
	free(head);
}

/* Returns elapsed wall-clock seconds since the previous call to periodt(). */
float periodt()
{
 clock_t marktime,temp;

	temp = oldmtime;
	marktime = clock();
	oldmtime = marktime;
	return	((marktime-temp)/CLK_TCK);
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

/* Builds an initial 2-partition greedily, picking each node from a
   randomized candidate list, using a doubly-linked list of unplaced
   nodes (sindex) to find candidates in O(n) per pass. */
void greedypart(int costa[],int nn,int ma[],
		int mb[],linknode sindex[],nodez alist[])
{
int nnhalf = nn/2, candid,z,head,point,numcand;
linknode *ind;
nodez *lptr;
register int x,y;
int cand[MAX_CAND_LIST_SIZE],candcost[MAX_CAND_LIST_SIZE],lowmax,lowmaxnode,highmin,highminnode;

	for (x = 1,ind = &sindex[1];x <=nn;x++,ind++)
		{
		ind -> back = x - 1;
		ind -> next = x + 1;
		}
	sindex[nn].next = 0;
	head = 1;
	for(z=0;z< cand_list_size;z++)
		{
		cand[z] = random(nn)+1;
		candcost[z] = 0;
		}
	numcand = cand_list_size;
	for (x = 1;x <= nnhalf; x++)
		{
	/* Put node from candidate list into set A and update cost */
                if (numcand > cand_list_size)
			numcand = cand_list_size;
		candid = ma[x] = cand[random(numcand)];
		if(candid == head)
			head = sindex[candid].next;
		else
			sindex[sindex[candid].back].next = sindex[candid].next;
		if (sindex[candid].next)
			sindex[sindex[candid].next].back = sindex[candid].back;
		for (y=1,z=alist[candid].node,lptr=alist[candid].next;y<=z;y++)
			{
			costa[lptr ->node]++;
			lptr = lptr ->next;
			}
		numcand = 0;
		for (z = 0; z<cand_list_size; z++)
			candcost[z] = 9999;
		highmin = 9999;
		highminnode = 0;
	/* Form candidate list for set B */
		point = head;
		while (point)
			{
			if(costa[point] < highmin)
				{
				numcand++;
				cand[highminnode] = point;
				candcost[highminnode] = costa[point];
				highmin = -9999;
				for (z = 0;z< cand_list_size;z++)
					if (candcost[z] > highmin)
						{
						highmin = candcost[z];
						highminnode = z;
						}
				}
			point = sindex[point].next;
			}
	/* Put node from candidate list into set B */
		if (numcand > cand_list_size)
			numcand = cand_list_size;
		candid = mb[x] = cand[random(numcand)];
		if(candid == head)
			head = sindex[candid].next;
		else
			sindex[sindex[candid].back].next = sindex[candid].next;
		if(sindex[candid].next)
			sindex[sindex[candid].next].back = sindex[candid].back;
		for (y=1,z=alist[candid].node,lptr=alist[candid].next;y<=z;y++)
			{
			costa[lptr ->node]--;
			lptr = lptr ->next;
			}
		for (z = 0; z < cand_list_size; z++)
			candcost[z] = -9999;
		lowmax = -9999;
		lowmaxnode = 0;
		numcand = 0;
	/* Form candidate list for set A */
		point = head;
		while (point)
			{
			if(costa[point] > lowmax)
				{
				numcand++;
				cand[lowmaxnode] = point;
				candcost[lowmaxnode] = costa[point];
				lowmax = 9999;
				for (z = 0;z< cand_list_size;z++)
					if (candcost[z] < lowmax)
						{
						lowmax = candcost[z];
						lowmaxnode = z;
						}
				}
			point = sindex[point].next;
			}
		}
	 /*  GREEDY PART */
/*	 for(x = 1;x<=nnhalf;x++)
		printf("%d  %d\n",ma[x],mb[x]);
	 printf("\n\n"); */
}

/* The following for functions are used in maintaining the "heap o' gains"
   for the greedy partition. Heapa (the nodes for set A) are sorted
   largest on top and heapb smallest on top, the two pairs of functions are
   only slightly different */
void upaheap(heapslot heap[],int i,crossindex index[])
{
int tempn,tempa,parent;

	while (i > 1)
		{
		parent = i/2;
		if (heap[i].alpha > heap[parent].alpha)
			{
			index[heap[parent].hnode].apoint = i;
			index[heap[i].hnode].apoint = parent;
			tempn = heap[parent].hnode;
			tempa = heap[parent].alpha;
			heap[parent].hnode = heap[i].hnode;
			heap[parent].alpha = heap[i].alpha;
			heap[i].hnode = tempn;
			heap[i].alpha = tempa;
			i = parent;
			}
		else
			i = 1;
		}
}

void downaheap(heapslot heap[],int i,int nn,crossindex index[])
{
int child,tempn,tempa;

	while (i < nn)
		{
		child = i << 1;
		if (child+1 <=nn && (heap[child+1].alpha > heap[child].alpha))
			child++;
		if(child <= nn && heap[i].alpha < heap[child].alpha)
			{
			index[heap[child].hnode].apoint = i;
			index[heap[i].hnode].apoint = child;
			tempn = heap[child].hnode;
			tempa = heap[child].alpha;
			heap[child].hnode = heap[i].hnode;
			heap[child].alpha = heap[i].alpha;
			heap[i].hnode = tempn;
			heap[i].alpha = tempa;
			i = child;
			}
		else
			i = nn;
		}
}

void upbheap(heapslot heap[],int i,crossindex index[])
{
int tempn,tempa,parent;

	while (i > 1)
		{
		parent = i/2;
		if (heap[i].alpha < heap[parent].alpha)
			{
			index[heap[parent].hnode].bpoint = i;
			index[heap[i].hnode].bpoint = parent;
			tempn = heap[parent].hnode;
			tempa = heap[parent].alpha;
			heap[parent].hnode = heap[i].hnode;
			heap[parent].alpha = heap[i].alpha;
			heap[i].hnode = tempn;
			heap[i].alpha = tempa;
			i = parent;
			}
		else
			i = 1;
		}
}

void downbheap(heapslot heap[],int i,int nn,crossindex index[])
{
int child,tempn,tempa;

	while (i < nn)
		{
		child = i << 1;
		if (child+1 <=nn && (heap[child+1].alpha < heap[child].alpha))
			child++;
		if(child <= nn && heap[i].alpha > heap[child].alpha)
			{
			index[heap[child].hnode].bpoint = i;
			index[heap[i].hnode].bpoint = child;
			tempn = heap[child].hnode;
			tempa = heap[child].alpha;
			heap[child].hnode = heap[i].hnode;
			heap[child].alpha = heap[i].alpha;
			heap[i].hnode = tempn;
			heap[i].alpha = tempa;
			i = child;
			}
		else
			i = nn;
		}
}

/* Builds an initial 2-partition greedily, same as greedypart but using
   two gain heaps (heapa/heapb) instead of a linear scan to pick each
   next candidate. */
void heappart(int costa[],int nn,int ma[],int mb[],nodez alist[])
{
int nnhalf = nn/2,top,cand,maxa,c,temptop,temp,clist[MAX_CAND_LIST_SIZE+1],temptr;
int z,y,candp;
register int x,csize;
heapslot templist[30],*heapa,*heapb,*ha,*hb,*tp,*tp2;
crossindex *index,*ind;
nodez *lptr;

	heapa = (heapslot *) calloc (nn+1,sizeof(heapslot));
	heapb = (heapslot *) calloc (nn+1,sizeof(heapslot));
	index = (crossindex *) calloc (nn+1,sizeof(crossindex));
	if (heapa == NULL || heapb == NULL || index == NULL)
		{
		puts("Couldn't allocate array in heappart");
		exit(0);
		}
	for(x = 1,ha = &heapa[1],hb = &heapb[1],ind = &index[1];
				x <=nn; x++,ha++,hb++,ind++)
		{
		ha->hnode = hb->hnode = x;
		ind->apoint = ind->bpoint = x;
		}

	for (x = 1;x <= nnhalf; x++)
		{
		/* Find a candidate list from the top of the A heap */
		clist[0] = heapa[1].hnode;
		temptop = 1;
		csize = 0;
		tp = templist;
		for (c = 1;c < cand_list_size;c++)
			{
			if(temptop * 2 < nn)
				{
				ha = &heapa[temptop << 1];
				if(ha->alpha > -9000)
					{
					tp->hnode = ha->hnode;
					tp->alpha = ha->alpha;
					tp++;
					csize++;
					}
				ha++;
				if(ha->alpha > -9000)
					{
					tp->hnode = ha->hnode;
					tp->alpha = ha->alpha;
					tp++;
					csize++;
					}
				maxa = -9999;
				for(temp = 0,tp2 = templist;temp<csize;
								temp++,tp2++)
					if(tp2->alpha > maxa)
						{
						maxa = tp2->alpha;
						top = tp2->hnode;
						temptr = temp;
						}
				if (maxa == -9999)
					break;
				clist[c] = top;
				temptop = index[top].apoint;
				templist[temptr].alpha = -9999;
				}
			else
				break;
			}
		ma[x] = cand = clist[random(c)];
		candp = index[cand].apoint;
		heapa[candp].alpha = -9999;
		downaheap(heapa,candp,nn,index);
		candp = index[cand].bpoint;
		heapb[candp].alpha = 9999;
		downbheap(heapb,candp,nn,index);
	/* Put node from candidate list into set A and update cost */
		for (y=1,z=alist[cand].node,lptr=alist[cand].next;y<=z;y++)
			{
			costa[lptr ->node]++;
			heapa[index[lptr->node].apoint].alpha++;
			heapb[index[lptr->node].bpoint].alpha++;
			upaheap(heapa,index[lptr->node].apoint,index);
			downbheap(heapb,index[lptr->node].bpoint,nn,index);
			lptr = lptr ->next;
			}
	/* Do the b side */
		clist[0] = heapb[1].hnode;
		temptop = 1;
		csize = 0;
                tp = templist;
		for (c = 1;c < cand_list_size;c++)
			{
			if (temptop * 2 < nn)
				{
				hb = &heapb[temptop << 1];
				if(hb->alpha < 9000)
					{
					tp->hnode = hb->hnode;
					tp->alpha = hb->alpha;
					tp++;
					csize++;
					}
				hb++;
				if(hb->alpha < 9000)
					{
					tp->hnode = hb->hnode;
					tp->alpha = hb->alpha;
					tp++;
					csize++;
					}
				maxa = 9999;
				for(temp = 0,tp2 = templist;temp<csize;
								temp++,tp2++)
					if(tp2->alpha < maxa)
						{
						maxa = tp2->alpha;
						top = tp2->hnode;
						temptr = temp;
						}
				if (maxa == 9999)
					break;
				clist[c] = top;
				temptop = index[top].bpoint;
				templist[temptr].alpha = 9999;
				}
			else
				break;
			}
		mb[x] = cand = clist[random(c)];
		candp = index[cand].apoint;
		heapa[candp].alpha = -9999;
		downaheap(heapa,candp,nn,index);
		candp = index[cand].bpoint;
		heapb[candp].alpha = 9999;
		downbheap(heapb,candp,nn,index);
	/* Put node from candidate list into set A and update cost */
		for (y=1,z=alist[cand].node,lptr=alist[cand].next;y<=z;y++)
			{
			costa[lptr->node]--;
			heapa[index[lptr->node].apoint].alpha--;
			heapb[index[lptr->node].bpoint].alpha--;
			downaheap(heapa,index[lptr->node].apoint,nn,index);
			upbheap(heapb,index[lptr->node].bpoint,index);
			lptr = lptr->next;
			}
		}
		free (heapa);
		free (heapb);
		free (index);
	 /*  HEAP GREEDY PART */
/*	for(x = 1;x<=nnhalf;x++)
		printf("%d  %d\n",ma[x],mb[x]);
	 printf("\n\n"); */
}

/* Find out if two nodes are adjacent; assumes input was ordered */
int finda_weight(nodez alist[],int node1,int node2)
{
register int x,z;
nodez *lptr;

	for(x=1,z=alist[node1].node,lptr = alist[node1].next;
	    x<=z && lptr->node <= node2;x++)
		{
		if(lptr->node == node2)
			return(1);
		lptr = lptr->next;
		}
	return(0);
}

/* Find a slight (gain < 3) positive swap and then perform it, using
   adjacency list */
void aslightswap(int ma[],int mb[],int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,worst,temp,found,gain,worstswap;
int worstx, worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[ma[x]] < 0)
			{
			y = 1;
			while (y <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					-(finda_weight(alist,ma[x],mb[y])<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain < 3)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				y++;
				}
			}
		x++;
		}
	    y = 1;
	    while (y <= nnhalf && !worst)
		{
		if (costa[mb[y]] > 0)
			{
			x = 1;
			while (x <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					-(finda_weight(alist,ma[x],mb[y])<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain <3)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				x++;
				}
			}
		y++;
		}
	    if(found)
		{
		temp = ma[worstx];
		ma[worstx] = mb[worsty];
		mb[worsty] = temp;
		for(x=1,z=alist[ma[worstx]].node,lptr=alist[ma[worstx]].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[mb[worsty]].node,lptr=alist[mb[worsty]].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1;x <= nnhalf; x++)
	for (y = 1;y<= nnhalf; y++)
		*cval += finda_weight(alist,ma[x],mb[y]);
	/* printf("Cross value after SLIGHT swap is %d\n",*cval); */
}



/* New improved generic two-exchange, with more postprocessing power!
   This is a "first swap" postprocessor with some improvements -- it looks
   first at those exchanges that are most likely to have a postive gain.
   The idea is to avoid having to look throught the whole list (n^2) */
void hswap(int igraph[],int ma[],int mb[],
	   int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,temp,found,gain,addx;
register int x,y,z;
nodez *lptr;

	do
	    {
	    found = 0;
	    x = 1;
	    while (!found && (x <= nnhalf))
		{
		if (costa[ma[x]] < 0)
			{
			y = 1;
			while (!found && (y <= nnhalf))
				{
				gain =  costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)]) << 1);
				if (gain > 0)
					found++;
				else
					{
					y++;
					}
				}
			if (!found)
				{
				x++;
				}
			}
		else
			{
			x++;
			}
		}
	    if (!found)
		{
		y = 1;
		}
	    while (!found && (y <= nnhalf))
		{
		if (costa[mb[y]] > 0)
			{
			x = 1;
			while (!found && (x <= nnhalf))
				{
				gain =  costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)]) << 1);
				if (gain > 0)
					found++;
				else
					{
					x++;
					}
				}
			if (!found)
				{
				y++;
				}
			}
		else
			{
			y++;
			}
		}
	    if(found)
		{
		/* x gets reused as the loop counter below, so remember
		   the found position in set A before that happens. */
		addx = x;
		temp = ma[addx];
		ma[addx] = mb[y];
		mb[y] = temp;
                for(x=1,z=alist[ma[addx]].node,lptr = alist[ma[addx]].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[mb[y]].node,lptr=alist[mb[y]].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
	/* printf("Cross value after HAL swap is %d\n",*cval); */
    *cval = 0;
    for (x = 1;x <= nnhalf; x++)
	for (y = 1;y<= nnhalf; y++)
		*cval += igraph[ind(ma[x],mb[y],nn)];

}

/* Find a slight (gain < 3) positive swap and then perform it */
void slightswap(int igraph[],int ma[],int mb[],int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,worst,temp,found,gain,worstswap;
int worstx, worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[ma[x]] < 0)
			{
			y = 1;
			while (y <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)])<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain < 3)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				y++;
				}
			}
		x++;
		}
	    y = 1;
	    while (y <= nnhalf && !worst)
		{
		if (costa[mb[y]] > 0)
			{
			x = 1;
			while (x <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)])<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain <3)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				x++;
				}
			}
		y++;
		}
	    if(found)
		{
		temp = ma[worstx];
		ma[worstx] = mb[worsty];
		mb[worsty] = temp;
		for(x=1,z=alist[ma[worstx]].node,lptr=alist[ma[worstx]].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[mb[worsty]].node,lptr=alist[mb[worsty]].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1;x <= nnhalf; x++)
	for (y = 1;y<= nnhalf; y++)
		*cval += igraph[ind(ma[x],mb[y],nn)];
	/* printf("Cross value after SLIGHT swap is %d\n",*cval); */
}

/* Find the slightest (gain == 1) positive swap and then perform it */
void slightestswap(int igraph[],int ma[],int mb[],
		   int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,worst,temp,found,gain,worstswap;
int worstx, worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[ma[x]] < 0)
			{
			y = 1;
			while (y <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)])<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain == 1)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				y++;
				}
			}
		x++;
		}
	    y = 1;
	    while (y <= nnhalf && !worst)
		{
		if (costa[mb[y]] > 0)
			{
			x = 1;
			while (x <= nnhalf && !worst)
				{
				gain = costa[mb[y]] - costa[ma[x]]
					- ((igraph[ind(ma[x],mb[y],nn)])<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain == 1)
						worst++;
					found++;
					worstswap = gain;
					worstx = x;
					worsty = y;
					}
				x++;
				}
			}
		y++;
		}
	    if(found)
		{
		temp = ma[worstx];
		ma[worstx] = mb[worsty];
		mb[worsty] = temp;
		for(x=1,z=alist[ma[worstx]].node,lptr=alist[ma[worstx]].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[mb[worsty]].node,lptr=alist[mb[worsty]].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1;x <= nnhalf; x++)
	for (y = 1;y<= nnhalf; y++)
		*cval += igraph[ind(ma[x],mb[y],nn)];
	/* printf("Cross value after SLIGHTEST swap is %d\n",*cval); */
}

