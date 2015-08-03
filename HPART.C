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
#include <time.h>
#include <string.h>

#define ind(i,j,nn) (i*nn+j)
#undef random
#define random(num) (rand() % (num))
#define hnode back
#define alpha next
#define apoint back
#define bpoint next

typedef struct anode{
	int node;
	struct anode *next;
	}nodez;

typedef nodez *nodep;

typedef struct {
	int back;
	int next;
	} indextype;

/* Function prototypes  */
void getgraph(int argc,char *argv[],
	      int  **igraph,int **a, int **b,indextype **sindex,
	      int *nn,int *ne,nodez **alist);
void readgraph(int ne,int nn,int igraph[],nodez alist[]);
void greedypart(int costa[],int nn,int ma[],
		int mb[],indextype sindex[],nodez alist[]);
void upaheap(indextype heap[],int i,indextype index[]);
void upbheap(indextype heap[],int i,indextype index[]);
void downaheap(indextype heap[],int i,int nn,indextype index[]);
void downbheap(indextype heap[],int i,int nn,indextype index[]);
void heappart(int costa[],int nn,int ma[],int mb[],nodez alist[]);
float periodt(void);
void remem(int **pa,int nn);
int finda_weight(nodez alist[],int node1,int node2);
void aslightswap(int ma[],int mb[],int costa[],int nn,int *cval,nodez alist[]);
void slightswap(int igraph[],int ma[],int mb[],int costa[],
		int nn,int *cval,nodez alist[]);
void slightestswap(int igraph[],int ma[],int mb[],int costa[],
		   int nn,int *cval,nodez alist[]);
void hswap(int igraph[],int ma[],int mb[],int costa[],int nn,
	   int *cval,nodez alist[]);

static FILE *f_in;
clock_t oldmtime;
int cand_list_size,big_flag = 1;

int main (int argc,char *argv[])
{
double sigx = 0.0 ,sigx2 = 0.0, sig;
int nn,ne,cval = 0,num_attemp = 0,mincval = 32600;
int *igraph, *ma,*mb,*costa,oneflag = 1;  /* ,*check,x; */
nodez *alist;
indextype *sindex;
clock_t startt;
float matcht = 0,partt= 0,swapt = 0,bestt,ttime,curtime,run_time;

	periodt();
	startt = clock();
	srand(32063);
	getgraph(argc,argv,&igraph,&ma,&mb,&sindex,
		 &nn,&ne,&alist);
	cand_list_size = atoi(argv[4]);
	if(argv[3][0] == '4')
		big_flag = 0;
	run_time = atof(argv[5]);

	readgraph(ne,nn,igraph,alist);
	curtime = clock()/CLK_TCK;
	while (curtime < run_time * 2)
		{
		num_attemp++;
		remem(&costa,nn);
		switch (argv[2][0])
			{
			case '1': heappart(costa,nn,ma,mb,alist);
				  partt += periodt();
				  break;
			case '2': greedypart(costa,nn,ma,mb,sindex,alist);
				  partt += periodt();
				  break;
			default : exit(0);
			}
		switch (argv[3][0])
			{
			case '1': hswap(igraph,ma,mb,costa,nn,&cval,alist);
				  break;
			case '2': slightswap(igraph,ma,mb,costa,nn,&cval,alist);
				  break;
			case '3': slightestswap(igraph,ma,mb,costa,nn,&cval,alist);
				  break;
			case '4': aslightswap(ma,mb,costa,nn,&cval,alist);
				  break;
			default : hswap(igraph,ma,mb,costa,nn,&cval,alist);
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
		/*	check = (int *) calloc (nn + 1,sizeof(int));
			for(x = 1; x<= nn / 2; x++)
				{
				check[ma[x]]++;
				check[mb[x]]++;
				}
			for (x = 1; x<= nn; x++)
				printf("\ncheck %d:  %d",x,check[x]);
			printf("\n\n");
			free(check);    */
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
//	printf("%f\n",ttime);
//	printf("%f\n",bestt);
//	printf("%f\n",matcht);
//	printf("%f\n",partt);
//	printf("%f\n",swapt);
//	printf("%lg\n",sigx);
//	printf("%lg\n",sigx2);
	printf("number of attempts = %d\n",num_attemp);
	free (alist);
	if (big_flag) free(igraph);
	free (ma); free(mb);
	free (sindex);
	return EXIT_SUCCESS;
}

/*  Inputs nn, ne, and allocates memory for matching.      */
void getgraph(int argc,char *argv[],
	      int **igraph,int **ma, int **mb, indextype **sindex,
	      int *nn,int *ne,nodez **alist)
{
char inputfile[12];
int nsqar,nn1,nnhalf;

	if (argc < 5)
		{
		puts("HPART: generates graph partition using RAGH");
		puts("Usage: HPART <inputfile> <modea> <modeb> <cl-size> <run-time>");
		puts("modea : 1    - heap greedy unmatched 2-partition");
		puts("        2    - greedy unmatched partition");
		puts("modeb   1    - first swap = 2-exchange");
		puts("        2    - slight swap (gain < 3)");
		puts("        3    - slightest swap (gain == 1)");
		puts("        4    - compact slight swap");
		exit(0);
		}
	strcpy(inputfile,argv[1]);
	f_in = fopen(inputfile,"r");
	if (f_in == NULL)
		{
		puts("Input file not found");
		exit(0);
		}
	fscanf(f_in,"%d,%d,%d",nn,ne);
        nn1 = (*nn)+1;
	nnhalf = (*nn) / 2;
	nsqar = (nn1) * (nn1);
	if(big_flag)
		{
		*igraph = (int *) calloc (nsqar,sizeof(int));
		if (*igraph == NULL)
			{
			puts("Couldn't make nxn array");
			exit(0);
			}
		}
	*ma = (int *) calloc (nnhalf + 1,sizeof(int));
	*mb = (int *) calloc (nnhalf + 1,sizeof(int));
	*sindex = (indextype *) calloc (nn1,sizeof(indextype));
	*alist = (nodez *) calloc (nn1,sizeof(nodez));
	if (*ma == NULL || *mb == NULL || *sindex == NULL || *alist == NULL)
		{
		puts("Couldn't allocate an array");
		exit(0);
		}
}

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

float periodt()
{
 clock_t marktime,temp;

	temp = oldmtime;
	marktime = clock();
	oldmtime = marktime;
	return	((marktime-temp)/CLK_TCK);
}

void remem(int **pa,int nn)
{
	*pa = (int *) calloc (nn+1,sizeof(int));
	if (*pa == NULL)
		{
		puts("Couldn't allocate an array");
		exit(0);
		}
}

void greedypart(int costa[],int nn,int ma[],
		int mb[],indextype sindex[],nodez alist[])
{
int nnhalf = nn/2, candid,*i,*j,z,head,point,numcand;
indextype *ind;
nodez *lptr;
register int x,y;
int cand[9],candc[9],lowmax,lowmaxnode,highmin,highminnode;

	for (x = 0,ind = &sindex[1];x <=nn;ind++)
		{
		ind -> back = x;
		ind -> next = ++x + 1;
		}
	sindex[nn].next = 0;
	head = 1;
	for(z=0;z< cand_list_size;z++)
		{
		cand[z] = random(nn)+1;
		candc[z] = 0;
		}
	numcand = cand_list_size;
	for (x = 1,i = &ma[1],j = &mb[1];x <= nnhalf; x++,i++,j++)
		{
	/* Put node from candidate list into set A and update cost */
                if (numcand > cand_list_size)
			numcand = cand_list_size;
		candid = *i = cand[random(numcand)];
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
			candc[z] = 9999;
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
				candc[highminnode] = costa[point];
				highmin = -9999;
				for (z = 0;z< cand_list_size;z++)
					if (candc[z] > highmin)
						{
						highmin = candc[z];
						highminnode = z;
						}
				}
			point = sindex[point].next;
			}
	/* Put node from candidate list into set B */
		if (numcand > cand_list_size)
			numcand = cand_list_size;
		candid = *j = cand[random(numcand)];
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
			candc[z] = -9999;
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
				candc[lowmaxnode] = costa[point];
				lowmax = 9999;
				for (z = 0;z< cand_list_size;z++)
					if (candc[z] < lowmax)
						{
						lowmax = candc[z];
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
void upaheap(indextype heap[],int i,indextype index[])
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

void downaheap(indextype heap[],int i,int nn,indextype index[])
{
int child,tempn,tempa;

	while (i < nn)
		{
		child = i << 1;
		if (child+1 <=nn && (heap[child+1].alpha > heap[child].alpha))
			child++;
		if(heap[i].alpha < heap[child].alpha && child <= nn)
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

void upbheap(indextype heap[],int i,indextype index[])
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

void downbheap(indextype heap[],int i,int nn,indextype index[])
{
int child,tempn,tempa;

	while (i < nn)
		{
		child = i << 1;
		if (child+1 <=nn && (heap[child+1].alpha < heap[child].alpha))
			child++;
		if(heap[i].alpha > heap[child].alpha && child <= nn)
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

void heappart(int costa[],int nn,int ma[],int mb[],nodez alist[])
{
int nnhalf = nn/2,top,cand,maxa,c,temptop,*i,*j,temp,clist[10],temptr;
int z,y,candp;
register int x,csize;
indextype templist[30],*heapa,*heapb,*index,*ha,*hb,*ind,*tp,*tp2;
nodez *lptr;

	heapa = (indextype *) calloc (nn+1,sizeof(indextype));
	heapb = (indextype *) calloc (nn+1,sizeof(indextype));
	index = (indextype *) calloc (nn+1,sizeof(indextype));
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

	for (x = 1,i = &ma[1],j = &mb[1];x <= nnhalf; x++,i++,j++)
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
		*i = cand = clist[random(c)];
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
		*j = cand = clist[random(c)];
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
int *i, *j, *worstx, *worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;  i = &ma[1];
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[*i] < 0)
			{
			y = 1; j = &mb[1];
			while (y <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					-(finda_weight(alist,*i,*j)<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain < 3)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				y++; j++;
				}
			}
		x++; i++;
		}
	    y = 1; j = &mb[1];
	    while (y <= nnhalf && !worst)
		{
		if (costa[*j] > 0)
			{
			x = 1; i = &ma[1];
			while (x <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					-(finda_weight(alist,*i,*j)<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain <3)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				x++;
				i++;
				}
			}
		y++;
		j++;
		}
	    if(found)
		{
		temp = *worstx;
		*worstx = *worsty;
		*worsty = temp;
		for(x=1,z=alist[*worstx].node,lptr=alist[*worstx].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[*worsty].node,lptr=alist[*worsty].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1,i= &ma[1];x <= nnhalf; x++,i++)
	for (y = 1, j= &mb[1];y<= nnhalf; y++,j++)
		*cval += finda_weight(alist,*i,*j);
	/* printf("Cross value after SLIGHT swap is %d\n",*cval); */
}



/* New improved generic two-exchange, with more postprocessing power!
   This is a "first swap" postprocessor with some improvements -- it looks
   first at those exchanges that are most likely to have a postive gain.
   The idea is to avoid having to look throught the whole list (n^2) */
void hswap(int igraph[],int ma[],int mb[],
	   int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,temp,found,gain;
register int x,y,z;
int *i, *j;
nodez *lptr;

	do
	    {
	    found = 0;
	    x = 1; i = &ma[1];
	    while (!found && (x <= nnhalf))
		{
		if (costa[*i] < 0)
			{
			y = 1; j = &mb[1];
			while (!found && (y <= nnhalf))
				{
				gain =  costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)]) << 1);
				if (gain > 0)
					found++;
				else
					{
					y++;
					j++;
					}
				}
			if (!found)
				{
				x++;
				i++;
				}
			}
		else
			{
			x++;
			i++;
			}
		}
	    if (!found)
		{
		y = 1;
		j = &mb[1];
		}
	    while (!found && (y <= nnhalf))
		{
		if (costa[*j] > 0)
			{
			x = 1; i = &ma[1];
			while (!found && (x <= nnhalf))
				{
				gain =  costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)]) << 1);
				if (gain > 0)
					found++;
				else
					{
					x++;
					i++;
					}
				}
			if (!found)
				{
				y++;
				j++;
				}
			}
		else
			{
			y++;
			j++;
			}
		}
	    if(found)
		{
		temp = *i;
		*i = *j;
		*j = temp;
                for(x=1,z=alist[*i].node,lptr = alist[*i].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[*j].node,lptr=alist[*j].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
	/* printf("Cross value after HAL swap is %d\n",*cval); */
    *cval = 0;
    for (x = 1,i= &ma[1];x <= nnhalf; x++,i++)
	for (y = 1,j= &mb[1];y<= nnhalf; y++,j++)
		*cval += igraph[ind(*i,*j,nn)];

}

/* Find a slight (gain < 3) positive swap and then perform it */
void slightswap(int igraph[],int ma[],int mb[],int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,worst,temp,found,gain,worstswap;
int *i, *j, *worstx, *worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;  i = &ma[1];
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[*i] < 0)
			{
			y = 1; j = &mb[1];
			while (y <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)])<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain < 3)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				y++; j++;
				}
			}
		x++; i++;
		}
	    y = 1; j = &mb[1];
	    while (y <= nnhalf && !worst)
		{
		if (costa[*j] > 0)
			{
			x = 1; i = &ma[1];
			while (x <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)])<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain <3)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				x++;
				i++;
				}
			}
		y++;
		j++;
		}
	    if(found)
		{
		temp = *worstx;
		*worstx = *worsty;
		*worsty = temp;
		for(x=1,z=alist[*worstx].node,lptr=alist[*worstx].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[*worsty].node,lptr=alist[*worsty].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1,i= &ma[1];x <= nnhalf; x++,i++)
	for (y = 1, j= &mb[1];y<= nnhalf; y++,j++)
		*cval += igraph[ind(*i,*j,nn)];
	/* printf("Cross value after SLIGHT swap is %d\n",*cval); */
}

/* Find the slightest (gain == 1) positive swap and then perform it */
void slightestswap(int igraph[],int ma[],int mb[],
		   int costa[],int nn,int *cval,nodez alist[])
{
int nnhalf = nn / 2,worst,temp,found,gain,worstswap;
int *i, *j, *worstx, *worsty,z;
register int x,y;
nodez *lptr;

	do
	    {
	    x = 1;  i = &ma[1];
	    worstswap =5000;
	    found = 0;
	    worst = 0;
	    while (x <= nnhalf && !worst)
		{
		if (costa[*i] < 0)
			{
			y = 1; j = &mb[1];
			while (y <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)])<<1);
				if (gain > 0 && gain < worstswap)
					{
					if (gain == 1)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				y++; j++;
				}
			}
		x++; i++;
		}
	    y = 1; j = &mb[1];
	    while (y <= nnhalf && !worst)
		{
		if (costa[*j] > 0)
			{
			x = 1; i = &ma[1];
			while (x <= nnhalf && !worst)
				{
				gain = costa[*j] - costa[*i]
					- ((igraph[ind(*i,*j,nn)])<<1);
				if (gain > 0 && gain< worstswap)
					{
					if (gain == 1)
						worst++;
					found++;
					worstswap = gain;
					worstx = i;
					worsty = j;
					}
				x++;
				i++;
				}
			}
		y++;
		j++;
		}
	    if(found)
		{
		temp = *worstx;
		*worstx = *worsty;
		*worsty = temp;
		for(x=1,z=alist[*worstx].node,lptr=alist[*worstx].next;x<=z;x++)
			{
			costa[lptr ->node] += 2;
			lptr = lptr->next;
			}
		for(x=1,z=alist[*worsty].node,lptr=alist[*worsty].next;x<=z;x++)
			{
			costa[lptr ->node] -= 2;
			lptr = lptr->next;
			}
		}
	    } while (found);
    *cval = 0;
    for (x = 1,i= &ma[1];x <= nnhalf; x++,i++)
	for (y = 1, j= &mb[1];y<= nnhalf; y++,j++)
		*cval += igraph[ind(*i,*j,nn)];
	/* printf("Cross value after SLIGHTEST swap is %d\n",*cval); */
}

