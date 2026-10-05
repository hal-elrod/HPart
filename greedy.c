/* greedy.c: partition construction (greedypart/heappart) and the
   local-search swap strategies (aslightswap/hswap/slightswap/
   slightestswap) for hpart.c. */

#include <stdio.h>
#include <stdlib.h>
#include "hpart.h"

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
			candcost[z] = POS_INF;
		highmin = POS_INF;
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
				highmin = NEG_INF;
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
			candcost[z] = NEG_INF;
		lowmax = NEG_INF;
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
				lowmax = POS_INF;
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
   only slightly different. Only heappart() calls these. */
static void upaheap(heapslot heap[],int i,crossindex index[])
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

static void downaheap(heapslot heap[],int i,int nn,crossindex index[])
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

static void upbheap(heapslot heap[],int i,crossindex index[])
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

static void downbheap(heapslot heap[],int i,int nn,crossindex index[])
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
				if(ha->alpha > -REMOVED_SLOT_MARGIN)
					{
					tp->hnode = ha->hnode;
					tp->alpha = ha->alpha;
					tp++;
					csize++;
					}
				ha++;
				if(ha->alpha > -REMOVED_SLOT_MARGIN)
					{
					tp->hnode = ha->hnode;
					tp->alpha = ha->alpha;
					tp++;
					csize++;
					}
				maxa = NEG_INF;
				for(temp = 0,tp2 = templist;temp<csize;
								temp++,tp2++)
					if(tp2->alpha > maxa)
						{
						maxa = tp2->alpha;
						top = tp2->hnode;
						temptr = temp;
						}
				if (maxa == NEG_INF)
					break;
				clist[c] = top;
				temptop = index[top].apoint;
				templist[temptr].alpha = NEG_INF;
				}
			else
				break;
			}
		ma[x] = cand = clist[random(c)];
		candp = index[cand].apoint;
		heapa[candp].alpha = NEG_INF;
		downaheap(heapa,candp,nn,index);
		candp = index[cand].bpoint;
		heapb[candp].alpha = POS_INF;
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
				if(hb->alpha < REMOVED_SLOT_MARGIN)
					{
					tp->hnode = hb->hnode;
					tp->alpha = hb->alpha;
					tp++;
					csize++;
					}
				hb++;
				if(hb->alpha < REMOVED_SLOT_MARGIN)
					{
					tp->hnode = hb->hnode;
					tp->alpha = hb->alpha;
					tp++;
					csize++;
					}
				maxa = POS_INF;
				for(temp = 0,tp2 = templist;temp<csize;
								temp++,tp2++)
					if(tp2->alpha < maxa)
						{
						maxa = tp2->alpha;
						top = tp2->hnode;
						temptr = temp;
						}
				if (maxa == POS_INF)
					break;
				clist[c] = top;
				temptop = index[top].bpoint;
				templist[temptr].alpha = POS_INF;
				}
			else
				break;
			}
		mb[x] = cand = clist[random(c)];
		candp = index[cand].apoint;
		heapa[candp].alpha = NEG_INF;
		downaheap(heapa,candp,nn,index);
		candp = index[cand].bpoint;
		heapb[candp].alpha = POS_INF;
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

/* Find out if two nodes are adjacent; assumes input was ordered. Only
   aslightswap() calls this (the big_flag-off / mode 4 path that avoids
   the igraph[] adjacency matrix). */
static int finda_weight(nodez alist[],int node1,int node2)
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
	    worstswap = WORST_SWAP_SENTINEL;
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
	    worstswap = WORST_SWAP_SENTINEL;
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
	    worstswap = WORST_SWAP_SENTINEL;
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
