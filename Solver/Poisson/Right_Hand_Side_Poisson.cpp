/*******************************************
 * Author: Michail Georgiou
<<<<<<< HEAD
 *  Last Modified: Time-stamp: <2014-05-01 12:34:39 mike_georgiou>
=======
 *  Last Modified: Time-stamp: <2014-04-28 15:54:00 mike_georgiou>
>>>>>>> f21e7fe8f645d773ef1fa9df298631c3fb0aefe6
 *
 *
Right_Hand_Side_Poisson.cpp -- This program computes the rhs of the poisson equation
of my code. The scheme that  I am using is presented in the
lessani-papalexandris paper, more specifically eq.66. In the vertical direction,
where the scheme is non-uniform i am using second order accurate scheme
*
* Written on Friday, 25 April 2014.
********************************************/
#include "Right_Hand_Side_Poisson.h"
#include "RHS.h"

void Right_Hand_Side_Poisson(double* rhs, double*** velocity_x,
                             double*** velocity_y, double*** velocity_z,
                             double*** rho_new, double*** rho,double*** rho_old,
                             double dx, double* dy, double dz,
                             double dt,
                             int ldx, int ldy, int ldz)
{
  double density_derivative, derivatives[3], divergence;


  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

        //Calculation of the time derivative of the Density

        // density_derivative =(3.0*rho_new[k][j][i]-4.0*rho[k][j][i]
        //                     +1.0*rho_old[k][j][i])/(2.*dt);

        //First Order, Just for Checking
        density_derivative = (rho_new[k][j][i]-rho[k][j][i])/(dt);

        //Divergence X-Component
        derivatives[0]=Divergence_X(velocity_x,dx,
                                    i, j, k);


        derivatives[1]=Divergence_Y(velocity_y,dy,
                                    i,j,k);


				derivatives[2]=Divergence_Z(velocity_z,dz,
																		i, j, k);

        divergence=0.;
        for (int vi=0; vi<3; vi++)
          divergence+=derivatives[vi];


        int index= A(i,ldx,
                     j, ldy,
                     k,ldz);

				 double result=(density_derivative+divergence)/dt;

				rhs[index+1] = result;


      }
    }
  }
}
