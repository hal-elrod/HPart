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
#include <string.h>
#include <stdio.h>

FILE *p_out,*p2_out;
int ne = 0;
int n,numn,a,b,i,j,w,twofiles = 0;
float p,r;
char outfile[12],outfile2[12];
long seed;
void help_me(void);

main(argc,argv)
int argc;
char *argv[];
{
double rando(long *seed);

   if (argc < 2 || argc > 3)
      help_me();
   strcpy(outfile,argv[1]);
   if (!(p_out = fopen(outfile,"w")))
      help_me();
   if (argc == 3)
      {
      strcpy(outfile2,argv[2]);
      if (!(p2_out = fopen(outfile2,"w")))
         help_me();
      twofiles = 1;
      }

   printf("Enter a seed for the random number generator (a big one): \n");
   scanf("%f",&seed);
   srand(seed);
   printf("Enter number of nodes, probability of edge\n");
   scanf("%d %f",&n,&p);
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
}

void help_me(void)
   {
   printf("HMAKE: a program to generate random graphs (no zero edges)\n");
   printf("Usage: HMAKE outfilename {outfile2}\n");
   exit(0);
   }
