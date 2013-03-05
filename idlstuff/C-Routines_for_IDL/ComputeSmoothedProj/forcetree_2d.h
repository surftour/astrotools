#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

void   force_costevaluate_2d(void) ;
int    force_getcost_single_2d(void);
int    force_getcost_quadru_2d(void);
void   force_resetcost_2d(void);
void   force_setupnonrecursive_2d(int no);
void   force_treeallocate_2d(int maxnodes, int maxpart);  
int    force_treebuild_2d(void);
int    force_treebuild_single_2d(int startnode, int *typelist, int *creatednodes);
void   force_treeevaluate_2d(int target, double one_over_s_of_a);
void   force_treeevaluate_direct_2d(int target, double one_over_s_of_a);
void   force_treeevaluate_single_2d(int tree, int targetpart, double epsilon);
void   force_treeevaluate_single_BH_2d(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential_2d(int target);
void   force_treeevaluate_potential_single_2d(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential_single_BH_2d(int tree, int targetpart, double epsilon);
void   force_treefree_2d(void);
void   force_update_node_2d(int no, int flag);
void   force_update_node_recursive_2d(int no);
void   force_update_size_of_parent_node_2d(int no);


float  INLINE_FUNC ngb_periodic_2d(float x);
float  ngb_select_closest_2d(int k, int n, float *arr, int *ind);
void   ngb_treeallocate_2d(int npart);
void   ngb_treebuild_2d(void);
float  ngb_treefind_2d(float xyz[3], int desngb, float hguess, int parttype, int **ngblistback, float **r2listback);
int    ngb_treefind_pairs_2d(float xyz[3], float hsml, int **ngblistback, float **r2listback);
int    ngb_treefind_variable_2d(float xyz[3], float hguess, int parttype, int **ngblistback, float **r2listback);
void   ngb_treefree_2d(void);
void   ngb_treesearch_2d(int);
void   ngb_treesearch_pairs_2d(int);
void   ngb_update_nodes_2d(void);






