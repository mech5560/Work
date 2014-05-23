#ifndef INTERMEDIATE_VELOCITY_H
#define INTERMEDIATE_VELOCITY_H

void Velocity_Residual_X( double*** residual_x, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y, double*** flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
			  double time_total,
                          int ldx, int ldy, int ldz);

void Velocity_Residual_Y( double*** residual_y,double*** velocity_x,
                          double*** velocity_y,double*** velocity_z,
                          double*** flux_x,double*** flux_y, double***flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
			  double time_total,
                          int ldx, int ldy, int ldz);


void Velocity_Residual_Z( double*** residual_z, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y, double*** flux_z,
                          double*** temperature, double Reynolds, double source,
                          double dx, double* dy, double dz,
			  double time_total,
                          int ldx, int ldy, int ldz);


#endif
