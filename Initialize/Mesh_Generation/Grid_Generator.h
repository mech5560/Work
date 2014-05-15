#ifndef GRID_GENERATOR_H
#define GRID_GENERATOR_H

void Hyperbolic_Mesh(int ldx, int ldy, int ldz,
                     double*** dy);

double *** Matrix_Alloc_New(int ldz, int ldy, int ldx,
                            int z_left,int z_right,
                            int y_left,int y_right,
                            int x_left,int x_right);

void Free_Matrix_New(double ***pA,
										 int z_left, int y_left, int x_left);


void Output_Data(int ldx, int ldy, int ldz,
                 double*** dy, char* name);
#endif
