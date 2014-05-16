/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-16 14:16:47 mike_georgiou>
 *
 *
Velocity_Residual_X.cpp -- This function computes the velocity residuals in the
X-direction (eq. 28 lessani-papalexandris paper). Still, fourth order accuracy
is used for the uniform grid direction (X and Z) and second order accuracy for
the Y-direction.
*
* Written on Thursday, 17 April 2014.
********************************************/
#include"Residuals.h"
#include"Residuals-inl.h"

void Velocity_Residual_X( double*** residual_x, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y, double*** flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
                          int ldx, int ldy, int ldz)
{

  //Variable for the calculation of the convective terms
  double  convection;

  //Variables for the calculation of the viscous terms
  double viscous_components[3],  viscous_total;



  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

        //Computing the convection term of the residuals by calling the

        //Convection Term
        convection= Convection_Term(velocity_x,
                                    flux_x, flux_y,flux_z,
                                    dx, dy, dz,
                                    i, j, k);

        //Viscous Term XX
        viscous_components[0] = Viscous_Component_XX(velocity_x,velocity_y,
                                                     velocity_z,
                                                     temperature, Reynolds,
                                                     dx, dy, dz,
                                                     i, j, k);


        //Viscous Term XY
        viscous_components[1] = Viscous_Component_XY(velocity_x,velocity_y,
                                                     temperature, Reynolds,
                                                     dx, dy,
                                                     i, j, k);

        //Viscous Term XY
        viscous_components[2] = Viscous_Component_XZ(velocity_x,velocity_z,
                                                     temperature, Reynolds,
                                                     dx, dz,
                                                     i, j, k);

        viscous_total=0.;
        for(int index=0; index<3; index++)
          viscous_total += viscous_components[index];

        residual_x[k][j][i] = -convection + viscous_total + source;

      }
    }
  }


}
