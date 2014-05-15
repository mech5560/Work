#ifndef RHS_H
#define RHS_H


double Divergence_X(double*** velocity_x, double dx,
                    int i, int j, int k);

double Divergence_Y(double*** velocity_y, double* dy,
										int i, int j, int k);
double Divergence_Z(double*** velocity_z, double dz,
                    int i, int j, int k);


#endif
