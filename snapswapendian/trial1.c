#include <stdio.h>
#include <stdlib.h>
#include <string.h>



int main(int argc,char *argv[])
{
  int i;
  char ch;
  char fname[100], new_fname[100];
  FILE *fin, *fout;


  /* --------- Check that the file was specified ------------- */
  if(argc!=2)
    {
      fprintf(stderr,"\n\nCall with input filename.\n\nConvertEndian <filename>");
      exit(1);
    }
  strcpy(fname, argv[1]);
  strcat(new_fname,".converted");

  /* --------- Make sure we can open the file ---------------- */
  if(!(fin=fopen(fname,"r")))
    {
      fprintf(stderr,"Can't open file '%s'.\n",fname);
      exit(1);
    }

  
  /* --------- Read in the header  ------------- */
   for (i=0; i<256; i++)
    {
      ch = getc(fin);
      putchar(ch);
    }

   printf("xxx\n\n");


   ch = getc(fin);
   while(!feof(fin))
     {
	i++;
	ch = getc(fin);
     }

   printf("the number of characters is: %d\n",i);

   fclose(fin);

   return 0;
}


































































