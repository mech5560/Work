/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-15 11:39:01 mike_georgiou>
 *
 *
Viscous_Component_XY.cpp -- This function computes the viscous component of the
X-Momentum equation
*
* Written on Wednesday, 23 April 2014.
********************************************/
#include"Residuals-inl.h"

double Viscous_Component_XY(double*** velocity_x, double*** velocity_y,
                            double*** temperature,double Reynolds,
                            double dx, double* dy,
                            int i, int j, int k)
{
  double  derivative_yy[2], derivative_yx[2][4];
  double dy_total, viscosity[4], viscous_terms[2];

  //calculation of the (du/dy)
  // j-1/2
  dy_total= dy[j]+dy[j-1];
  derivative_yy[0] =Derivative(velocity_x[k][j][i],
                               velocity_x[k][j-1][i],
                               dy_total, 1);

  // j+1/2
  dy_total= dy[j]+dy[j+1];
  derivative_yy[1] =Derivative(velocity_x[k][j+1][i],
                               velocity_x[k][j][i],
                               dy_total,1);




  // Calculation of the (dv/dx)
  //j-1/2
  // j-1
  derivative_yx[0][0] = 9./(8.)*Derivative(velocity_y[k][j-1][i+1],
                                           velocity_y[k][j-1][i-1],
                                           dx,2);

  derivative_yx[0][1] = -1./(8.)*Derivative(velocity_y[k][j-1][i+3],
                                            velocity_y[k][j-1][i-3],
                                            dx,6);

  // j
  derivative_yx[0][2] = 9./(8.)*Derivative(velocity_y[k][j][i+1],
                                           velocity_y[k][j][i-1],
                                           dx,2);

  derivative_yx[0][3] = -1./(8.)*Derivative(velocity_y[k][j][i+3],
                                            velocity_y[k][j][i-3],
                                            dx,6);

  //j+1/2
  // j
  derivative_yx[1][0] =derivative_yx[0][2];
  derivative_yx[1][1] =derivative_yx[0][3];

  // j+1
  derivative_yx[1][2] = 9./(8.)*Derivative(velocity_y[k][j+1][i+1],
                                           velocity_y[k][j+1][i-1],
                                           dx,2);

  derivative_yx[1][3] = -1./(8.)*Derivative(velocity_y[k][j+1][i+3],
                                            velocity_y[k][j+1][i-3],
                                            dx,6);

  //computing the interpolated derivative d/dy(dv/dx)
  double sum[2];
  double total_derivative_x[2];

  for (int vi=0; vi<2; vi++)
    {
      sum[0] = derivative_yx[vi][0]+derivative_yx[vi][1];
      sum[1] = derivative_yx[vi][2]+derivative_yx[vi][3];


      total_derivative_x[vi]= Interpolation_Y(sum[1], dy[j],
                                              sum[0], dy[j]);
    }


  // Computing the Viscosity
  for (int vi=0; vi<2; vi++)
    {
      viscosity[vi]=
        Viscosity_Calculator(Interpolation_Y(temperature[k][j+vi][i],
                                             dy[j+1],
                                             temperature[k][j+vi-1][i],
                                             dy[j]));
    }

  for (int vi=0; vi<2; vi++)
    {
      viscous_terms[vi]=viscosity[vi]*(derivative_yy[vi]
                                       +total_derivative_x[vi]);
    }

  //computing the viscous component in the y-direction
  dy_total=2.*dy[j];
  double viscous_component= 1./(Reynolds*dy_total)*(viscous_terms[1]-
                                                    viscous_terms[0]);


  return viscous_component;
}
