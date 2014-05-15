/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-04-30 10:12:21 mike_georgiou>
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


const double pi = 4.*atan(1.), a = .98346;

void Hyperbolic_Mesh(int ldx, int ldy, int ldz,
                     double*** dy)
{

  double atanha = .5*(log((1+a)/(1-a)));
  double ksi;


  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

        ksi =  -1. + 2.*(j)/(ldy);
        dy[k][j][i] = (1./a)*tanh(ksi*atanha);

      }
    }
  }





}
