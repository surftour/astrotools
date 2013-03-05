#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


#ifdef PEDANTIC_MEMORY_HANDLER

#define MAXBLOCKS 5000

#ifdef PEDANTIC_MEMORY_CEILING
static size_t TotBytes;
static void *Base;
#endif

static unsigned long Nblocks;
static double Highmark;


static void **Table;
static size_t *BlockSize;

void mymalloc_init(void)
{
#ifdef PEDANTIC_MEMORY_CEILING
  size_t n;
#endif

  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));

#ifdef PEDANTIC_MEMORY_CEILING
  n = PEDANTIC_MEMORY_CEILING * 1024.0 * 1024.0;

  if(!(Base = malloc(n)))
    {
      printf("Failed to allocate memory for `Base' (%d Mbytes).\n", (int) PEDANTIC_MEMORY_CEILING);
      endrun(122);
    }

  TotBytes = FreeBytes = n;
#endif

  AllocatedBytes = 0;
  Nblocks = 0;
  Highmark = 0;
}



void *mymalloc(size_t n)
{
  if((n % 8) > 0)
    n = (n / 8 + 1) * 8;

  if(n < 8)
    n = 8;

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in mymalloc().\n", ThisTask);
      endrun(813);
    }

#ifdef PEDANTIC_MEMORY_CEILING
  if(n > FreeBytes)
    {
      printf("Task=%d: Not enough memory in mymalloc(n=%g MB).  FreeBytes=%g MB\n",
	     ThisTask, n / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;
#else
  Table[Nblocks] = malloc(n);
  if(!(Table[Nblocks]))
    {
      printf("failed to allocate %g MB of memory. (presently allocated=%g MB)\n",
	     n / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));
      endrun(18);
    }
#endif

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;

  Nblocks += 1;

  /*
     if(AllocatedBytes / (1024.0 * 1024.0) > Highmark)
     {
     Highmark = AllocatedBytes / (1024.0 * 1024.0);
     printf("Task=%d:   new highmark=%g MB\n",
     ThisTask, Highmark);
     fflush(stdout);
     }
   */

  return Table[Nblocks - 1];
}


void *mymalloc_msg(size_t n, char *message)
{
  if((n % 8) > 0)
    n = (n / 8 + 1) * 8;

  if(n < 8)
    n = 8;

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in mymalloc().\n", ThisTask);
      endrun(813);
    }

#ifdef PEDANTIC_MEMORY_CEILING
  if(n > FreeBytes)
    {
      printf
	("Task=%d: Not enough memory in for allocating (n=%g MB) for block='%s'.  FreeBytes=%g MB  AllocatedBytes=%g MB.\n",
	 ThisTask, n / (1024.0 * 1024.0), message, FreeBytes / (1024.0 * 1024.0),
	 AllocatedBytes / (1024.0 * 1024.0));
      endrun(812);
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;
#else
  Table[Nblocks] = malloc(n);
  if(!(Table[Nblocks]))
    {
      printf("failed to allocate %g MB of memory for block='%s'. (presently allocated=%g MB)\n",
	     n / (1024.0 * 1024.0), message, AllocatedBytes / (1024.0 * 1024.0));
      endrun(18);
    }
#endif

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;

  Nblocks += 1;

  return Table[Nblocks - 1];
}




void myfree(void *p)
{
  if(Nblocks == 0)
    endrun(76878);

  if(p != Table[Nblocks - 1])
    {
      printf("Task=%d: Wrong call of myfree() - not the last allocated block!\n", ThisTask);
      fflush(stdout);
      endrun(814);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
#ifdef PEDANTIC_MEMORY_CEILING
  FreeBytes += BlockSize[Nblocks];
#else
  free(p);
#endif
}

void myfree_msg(void *p, char *msg)
{
  if(Nblocks == 0)
    endrun(76878);

  if(p != Table[Nblocks - 1])
    {
      printf("Task=%d: Wrong call of myfree() - '%s' not the last allocated block!\n", ThisTask, msg);
      fflush(stdout);
      endrun(8141);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
#ifdef PEDANTIC_MEMORY_CEILING
  FreeBytes += BlockSize[Nblocks];
#else
  free(p);
#endif
}



void *myrealloc(void *p, size_t n)
{
  if((n % 8) > 0)
    n = (n / 8 + 1) * 8;

  if(n < 8)
    n = 8;

  if(Nblocks == 0)
    endrun(76879);

  if(p != Table[Nblocks - 1])
    {
      printf("Task=%d: Wrong call of myrealloc() - not the last allocated block!\n", ThisTask);
      fflush(stdout);
      endrun(815);
    }

  AllocatedBytes -= BlockSize[Nblocks - 1];
#ifdef PEDANTIC_MEMORY_CEILING
  FreeBytes += BlockSize[Nblocks - 1];
#endif


#ifdef PEDANTIC_MEMORY_CEILING
  if(n > FreeBytes)
    {
      printf("Task=%d: Not enough memory in myremalloc(n=%g MB). previous=%g FreeBytes=%g MB\n",
	     ThisTask, n / (1024.0 * 1024.0), BlockSize[Nblocks - 1] / (1024.0 * 1024.0),
	     FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }
  Table[Nblocks - 1] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;
#else
  Table[Nblocks - 1] = realloc(Table[Nblocks - 1], n);
  if(!(Table[Nblocks - 1]))
    {
      printf("failed to reallocate %g MB of memory. previous=%g FreeBytes=%g MB\n",
	     n / (1024.0 * 1024.0), BlockSize[Nblocks - 1] / (1024.0 * 1024.0),
	     FreeBytes / (1024.0 * 1024.0));
      endrun(18123);
    }
#endif

  AllocatedBytes += n;
  BlockSize[Nblocks - 1] = n;

  return Table[Nblocks - 1];
}



#else

void mymalloc_init(void)
{
  AllocatedBytes = 0;
}

void *mymalloc(size_t n)
{
  void *ptr;

  if(n < 8)
    n = 8;

  ptr = malloc(n);

  if(!(ptr))
    {
      printf("failed to allocate %g MB of memory.\n", n / (1024.0 * 1024.0));
      endrun(14);
    }

  return ptr;
}


void *mymalloc_msg(size_t n, char *message)
{
  void *ptr;

  if(n < 8)
    n = 8;

  ptr = malloc(n);

  if(!(ptr))
    {
      printf("failed to allocate %g MB of memory when trying to allocate '%s'.\n",
	     n / (1024.0 * 1024.0), message);
      endrun(14);
    }

  return ptr;
}



void *myrealloc(void *p, size_t n)
{
  void *ptr;

  ptr = realloc(p, n);

  if(!(ptr))
    {
      printf("failed to re-allocate %g MB of memory.\n", n / (1024.0 * 1024.0));
      endrun(15);
    }

  return ptr;
}

void myfree(void *p)
{
  free(p);
}

void myfree_msg(void *p, char *msg)
{
  free(p);
}



#endif
