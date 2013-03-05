

void ngb2d_treebuild(float **pospointer,int Npart,int MinMaxFlag,float XminOpt[2],float XmaxOpt[2]);

void ngb2d_treefree(void);
void ngb2d_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */

float ngb2d_treefind(float xyz[2], int desngb, float hguess,int **ngblistback, float **r2listback,float hmax,int *ngbfound);


float ngb2d_treetest(float xyz[2], int desngb, float hguess,int **ngblistback, float **r2listback);
