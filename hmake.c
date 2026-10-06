/* ......................................................................
   HMAKE:  Generates random graphs in specified formats

   by Hal C. Elrod      2/89        University of Texas
   rewritten 3/90

   Format:
      HMAKE <outfile> {outfile2}

   Asks for:
        seed = random number seed
      n = number of nodes (must be even for partition)
      p = probability of an edge

   Returns:  If one outfile is specified, MAKEGR produces one file
        in the following format.

         n,   ne        <- ne = number of edges
         node1,   node2,   weight(1)
           .     .        .
           .      .        .
         node1,   node,   weight(ne)

        If both outfiles specified, MAKEGR makes both the above
        file and one in the input format for the Kernighan & Lin
        partitioning program.

   Outfile2:   F A B C P Q             <- flags
         M 1 C                   <- modules or nodes
           .
           .
         M n C
         L x y w         <- links, a<w<b
         E         <- end

.................................................................*/

#include <stdlib.h>
#include <stdio.h>
#include "hmake.h"

FILE *p_out,*p2_out;
int ne = 0;
int n,numn,i,j,twofiles = 0;
float p,r;
long seed;

/* Reads a seed, node count, and edge probability from stdin, then writes
   a random 0-1 graph (every edge weight 1) to argv[1] in "n,ne" + edge-list
   format, and optionally a second file in Kernighan-Lin input format
   to argv[2]. See the file header above for both formats. */
int main(int argc,char *argv[])
{
   if (argc < 2 || argc > 3)
      help_me();
   if (!(p_out = fopen(argv[1],"w")))
      help_me();
   if (argc == 3)
      {
      if (!(p2_out = fopen(argv[2],"w")))
         {
         fclose(p_out);
         help_me();
         }
      twofiles = 1;
      }

   printf("Enter a seed for the random number generator (a big one): \n");
   if (scanf("%ld",&seed) != 1)
      {
      printf("Error: couldn't read a seed value.\n");
      exit(1);
      }
   srand(seed);
   printf("Enter number of nodes, probability of edge\n");
   if (scanf("%d %f",&n,&p) != 2)
      {
      printf("Error: couldn't read node count / edge probability.\n");
      exit(1);
      }
   if (twofiles)
      {
      fputs("F A B C P Q\n",p2_out);
      for (numn= 1;numn <= n;numn++)
         fprintf(p2_out,"M %d C\n",numn);
      }
   fprintf(p_out,"                                   \n");
   for (i= 1;i<=n;i++)
      {
      for (j=i+1;j<=n;j++)
         {
         r = (float)rand()/(float)RAND_MAX;
         if(r<p)
         {
         ne++;
         fprintf(p_out,"%d,%d,1\n",i,j);
         if (twofiles )
            fprintf(p2_out,"L %d %d 1\n",i,j);
         }
         }
      }
   fseek(p_out,0,SEEK_SET);
   fprintf(p_out,"%d,%d",n,ne);
   fclose(p_out);
   if (twofiles)
      {
                fputs("E\n",p2_out);
      fclose(p2_out);
      }
   return 0;
}

/* Prints usage and exits; called on any bad argument or unopenable file. */
void help_me(void)
   {
   printf("HMAKE: a program to generate random graphs (no zero edges)\n");
   printf("Usage: HMAKE outfilename {outfile2}\n");
   exit(0);
   }
