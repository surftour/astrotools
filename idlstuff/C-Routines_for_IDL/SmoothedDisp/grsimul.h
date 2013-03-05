
void  calculate_force_on_x(double xpos[][3],double accel[][3],double pot[],int nchips);
int   read_neighbor_list(int board);
int   get_neighbor_list(int ichip,int *list);
void  initiliaze_grape3(void);
void  free_grape3(void);
void  set_scales(double mi,double ma);
void  set_mass_correction_factor(double c);
int   number_of_available_chips(void);
int   number_of_chips_per_board(void);
int   number_of_boards(void);
void  set_n(int n);
void  set_xj(int j, double x[3]);
void  set_mj(int j,double m);
void  set_eps2(double eps2);
void  set_eps2_to_chip(double eps2,int ichip);
void  set_h2(double h2,int chip);
int   max_particle_number(void);
int   max_neighbor_number(void);
