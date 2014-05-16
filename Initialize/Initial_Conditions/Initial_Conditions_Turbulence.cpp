/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-06 15:20:59 mike_georgiou>
 *
 *
Initialize_Conditions_Turbulence.cpp -- This function initializes the problem for the
turbulent pipe flow. Specifically, the streamwise velocity "w" is
being initialized by using the log law profile. Also, in all
directions, pertubations of order ??? will be introduced.

The pressure gradient will be kept constant in the streamwise direction.
*
* Written on Monday,  5 May 2014.
********************************************/
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <ctime>

#include"Initial_Conditions_Turbulence.h"

void Initial_Conditions_Turbulence(double*** velocity_x,
                                   double*** velocity_y,
                                   double*** velocity_z,
                                   double Reynolds, double dt,
                                   double dx, double* dy, double dz,
                                   int ldx, int ldy, int ldz)
{
  // Reading data for the initialization of the turbulent flow
  double Reynolds_wall=0.;
  Turbulence_Reader(&Reynolds_wall);

  // Computing the friction velocity and the viscosity which is
  // necessary for the implementation of the law of the wall
  double friction_velocity=1.;
  double viscosity = friction_velocity/Reynolds_wall;


  // Implementing the law of the wall
  for (int k=0; k<ldz; k++){
    double y_local=0.;
    for (int j=0; j<ldy/2 +1; j++){

      //defining the vertical position that we are
      y_local += 2.*dy[j];

      for (int i=0; i<ldx; i++){
	double condition=y_local*friction_velocity/viscosity;

        if(condition<10.)
          {
            velocity_z[k][j][i]=y_local*friction_velocity/viscosity;
          }
        else
          {
            velocity_z[k][j][i]=2.5*log (y_local*friction_velocity/viscosity)+5.;
          }

        //Reflecting the initial condition
        velocity_z[k][ldy-1-j][i]= velocity_z[k][j][i];

      }
    }
  }


  //Seeding the rand() function
  srand(time(NULL));


  double reference_velocity=velocity_z[0][ldy/2+1][0];
  double rnd;
  //Introducing pertubations
  for (int k = 0; k < ldz; k++){
    for (int j = 0; j < ldy; j++){
      for (int i = 0; i < ldx; i++){

  	//rnd=rand()%100000 -50000.;
        velocity_x[k][j][i]= ( (velocity_x[k][j][i] +
  				1e-7*reference_velocity*(rand()%100000 -50000.))/
  			       reference_velocity);

  	//rnd=rand()%100000 -50000.;
        velocity_y[k][j][i]= ( (velocity_y[k][j][i] +
  				1e-7*reference_velocity*(rand()%100000 -50000.))/
  			       reference_velocity);

        //rnd=rand()%100000 -50000.;
        velocity_z[k][j][i]=( (velocity_z[k][j][i] +
  			       1e-7*reference_velocity*(rand()%100000 -50000.))/
  			      reference_velocity);
      }
    }
  }


}
