#ifndef INITIALIZE_ALLOCATE_H
#define INITIALIZE_ALLOCATE_H
double *** Matrix_Allocator(int ldz, int ldy, int ldx,
                            int z_left,int z_right,
                            int y_left,int y_right,
                            int x_left,int x_right);

void Free_Matrix(double ***pA,
								 int z_left, int y_left, int x_left);




struct Arrays{
  double
    *dy,

    ***velocity_x, ***velocity_y, ***velocity_z,
    ***velocity_x_new, ***velocity_y_new, ***velocity_z_new,
    ***velocity_x_tilda, ***velocity_y_tilda, ***velocity_z_tilda,

    ***flux_x, ***flux_y, ***flux_z,
    ***flux_x_new, ***flux_y_new, ***flux_z_new,

    ***rho, ***rho_new, ***rho_old,

    ***pressure, 

    ***temperature, ***temperature_new, ***temperature_avg,

    ***residual_y, ***residual_z, ***residual_x,
    ***residual_x_old, ***residual_y_old, ***residual_z_old;
};

typedef   struct Arrays Ar;



#endif

