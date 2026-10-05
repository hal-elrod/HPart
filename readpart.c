/* readpart.c: reads the input graph file for hpart.c.

   Input file format:
	number of nodes, number of edges
	node, node, weight
	 .      .     .
	 .	.     .
	node, node, weight
   ....................................................................... */

#include <stdio.h>
#include <stdlib.h>
#include "hpart.h"

static FILE *f_in;

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
