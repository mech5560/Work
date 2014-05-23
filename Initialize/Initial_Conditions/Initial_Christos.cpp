/*******************************************
 * Author: Michail Georgiou 
*  Last Modified: Time-stamp: <2014-05-22 14:55:38 mike_georgiou>   
*
*
Initial_Christos.cpp -- This function introduces as initial
conditions the expected solution of the 2D problem that Christos gave
me. The resonable result is that the poisson solver will reach
converge very fast.
*
* Written on Wednesday, 21 May 2014.
********************************************/
#include <math.h>
#include <iostream>
#include <stdio.h>
using namespace std;


void Initial_Christos(double*** velocity_x, double*** velocity_y,
		      double*** velocity_z,
		      double dx,
		      int ldx, int ldy, int ldz)
{
  //  double omega = 1 +sin(2.*pi*time*time);
  static const double pi=4.*atan(1.);


for (int k = 0; k < ldz; k++){
  double y_local =0.;
  for (int j = 0; j < ldy; j++){

    y_local+=dx/2.;
    //Definition of the y components of the forcing term
    double y_component_1 = 3.*y_local*y_local-2.*y_local;


    double y_component_2 = y_local*y_local*(y_local-1.);

    double x_local=0.;
    for (int i = 0; i < ldx; i++){
      
      x_local += dx/2.;


      velocity_x[k][j][i]= cos(2.*pi*(x_local))*y_component_1;


      velocity_y[k][j][i] =
		2.*pi*sin(2.*pi*(x_local))*y_component_2;
      
      velocity_z[k][j][i] =0.;
      x_local += dx/2.;
      
    }
    y_local+=dx/2.;
  }	    
 }   



}
