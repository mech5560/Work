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




inline double Derivative(double pressure_right, double pressure_left,
                         double dx, int order)
{

  return 1./(order * dx)*( (pressure_right) - (pressure_left));

}


inline double Interpolation_Y(double value_right, double dy_right,
                              double value_left, double dy_left)
{
  double coefficient[2];
  coefficient[0] = dy_left/(dy_left+dy_right); // Right coeff
  coefficient[1] = dy_right/(dy_left+dy_right); // left coeff

  return coefficient[0]*value_right+
    coefficient[1]*value_left;
}





#endif

