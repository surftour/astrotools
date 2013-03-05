#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "globvars.h"


int main(int argc, char **argv)
{
  FILE *fd;
  char input_fname[200];
  int k, dummy, n_in_blk;
  int dummy_total=0;
  int fnum=0;
  float junk;
  char *p;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  if(argc<2)
    {
      printf("\nReadIC <IC/snapshot filename> \n\n");
      exit(0);
    }

  strcpy(input_fname, argv[1]);


	/* =============================================== */

        if(!(fd=fopen(input_fname,"r")))
        {
          printf("can't open file `%s`\n",input_fname);
          exit(0);
        }

	SKIP;
	printf("%10d\tHeader\n",dummy);  dummy_total+= dummy+8;
	fread(&header1, sizeof(header1), 1, fd);
	SKIP;

	/* while (!feof(fd)) */
	do
	{
	   SKIP;
	   n_in_blk=dummy/4;
	   dummy_total+= dummy+8;
	   printf("%10d\tField #%2d   (n=%6d)\n",dummy,fnum,n_in_blk);
	   for(k=0;k<n_in_blk;k++) {
		fread(&junk, sizeof(float), 1, fd);
		/* if ((k % 100000) == 0) printf("%d ... ",k); */
		}
	   /* printf("\n"); */
	   SKIP;
	   fnum++;
	} while (!feof(fd));

	fclose(fd);

	printf("----------\n");
	printf("%10d\tTotal Bytes (included padding)\n\n",dummy_total);


	/* =============================================== */


  return(0);
}






