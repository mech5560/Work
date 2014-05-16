/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-16 09:39:22 mike_georgiou>
 *
 *
Hyperbolic_Mesh.cpp -- This function generates non-uniform coordinates for the
vertical direction. To generate these coordinates I used information from the
code of Bamdad Lessani.

The input of this code are the number of solution points in the three directions
and the expansion of the number of ghost cells in all the edges of the domain.

The output of this code is a matrix with the y coordinates of my domain.

*
* Written on Wednesday 30 April 2014.
********************************************/

#include <cmath>
#include <iostream>
using namespace std;

const double pi = 4.*atan(1.), a = .98346;

void Hyperbolic_Mesh(double* dy, double *y,
		     double length_y,
		     int ldy, int ly, int ry)
{

  //  double dy_temp[ldy+1];
  double atanha = .5*(log((1+a)/(1-a)));
  double ksi;


  for (int j=0; j<=ldy; j++)
    {

      ksi =  -1. + 2.*(j)/(ldy);
      y[j] = (1./a)*tanh(ksi*atanha);
 
    }
   
  for (int j=0; j<ldy; j++)
    {
      dy[j] = 0.5*(y[j+1]-y[j]);
    }

  double y_temp = -1.;
  for (int j=0; j<ldy; j++)
    {
      y_temp += dy[j];
      y[j]=y_temp;
      y_temp += dy[j];
    }

	//Boundary Conditions			
	dy[-1] = dy[0];
	dy[ldy] = dy[ldy-1];

	y[-1] = y[0];
	y[ldy] = y[ldy-1];

}
