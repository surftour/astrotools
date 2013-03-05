#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct ray_source_struct {
	int ray_source_ID;
	float *ray_nh_list;
	float *ray_Z_list;
} ray_sources;


int main()
{
	char out_file_2[100];
		printf("\nFilename for integrated values (NH, etc): ");
		scanf("%s",out_file_2);	
		FILE *output_file;
    	if ((output_file = fopen(out_file_2, "rb")) == NULL) 
    		fprintf(stderr, "Cannot open %s\n", "output_file");	
	
	float simulation_time;
	int ray_origin_type;
	int number_of_ray_origins;
	int number_of_rays_per_origin;

    fread(&simulation_time,sizeof(float),1,output_file);
    fread(&ray_origin_type,sizeof(int),1,output_file);
    fread(&number_of_ray_origins,sizeof(int),1,output_file);
    fread(&number_of_rays_per_origin,sizeof(int),1,output_file);
	printf("%i %i %i %f \n",ray_origin_type,number_of_ray_origins,number_of_rays_per_origin,simulation_time);
    
    int nO = number_of_ray_origins;
    int nR = number_of_rays_per_origin;
    int i=0,j=0;
    float *theta_list;
    theta_list = malloc(nR * sizeof(float));
    float *phi_list;
    phi_list = malloc(nR * sizeof(float));
	for (i=0; i<nR; i++) {
    	fread(&theta_list[i],sizeof(float),1,output_file);
   		fread(&phi_list[i],sizeof(float),1,output_file);
	}
	
	ray_sources *ray_origins;
	ray_origins = malloc(nO*sizeof(ray_sources));
	for (i=0; i<nO; i++) {
	    fread(&ray_origins[i].ray_source_ID,sizeof(int),1,output_file);
		ray_origins[i].ray_nh_list = malloc(nR * sizeof(float));
		ray_origins[i].ray_Z_list = malloc(nR * sizeof(float));
		for (j=0; j<nR; j++) {
//		for (j=0; j<10; j++) {
	    	fread(&ray_origins[i].ray_nh_list[j],sizeof(float),1,output_file);
	    	fread(&ray_origins[i].ray_Z_list[j],sizeof(float),1,output_file);
			printf("%e %e \n",ray_origins[i].ray_nh_list[j],ray_origins[i].ray_Z_list[j]);
		}
	}

    fclose(output_file);
	return 0;	
}
