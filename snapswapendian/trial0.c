#include <stdio.h>
#include <stdlib.h>



int main()
{

      char c;

      printf("float %3d\n", sizeof(float));
      printf("double %3d\n", sizeof(double));
      printf("long double %3d\n", sizeof(long double));
      printf("short int %3d\n", sizeof(short int));
      printf("int %3d\n", sizeof(int));
      printf("long %3d\n", sizeof(long));
      printf("string %3d\n", sizeof("string"));
      printf("char %3d\n", sizeof(c));

      return(0);

}


