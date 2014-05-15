/*******************************************
 * Author: Michail Georgiou
<<<<<<< HEAD
*  Last Modified: Time-stamp: <2014-05-15 12:48:16 mike_georgiou>
=======
*  Last Modified: Time-stamp: <2014-04-28 16:00:10 mike_georgiou>
>>>>>>> f21e7fe8f645d773ef1fa9df298631c3fb0aefe6
*
*
Velocity_Residual_Z.cpp -- This function computes the velocity residuals in the
Z-direction (eq. 28 lessani-papalexandris paper). Still, fourth order accuracy
is used for the uniform grid direction (X and Z) and second order accuracy for
the Y-direction.
*
* Written on Thursday, 17 April 2014.
********************************************/


#include"../../Header_Files/Data.h"

#include"Residuals.h"
#include"Residuals-inl.h"

void Velocity_Residual_Z( double*** residual_z, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y, double*** flux_z,
                          double*** temperature, double Reynolds, double source,
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


        //Convection Term
        convection= Convection_Term(velocity_z,
                                    flux_x, flux_y,flux_z,
                                    dx, dy, dz,
                                    i, j, k);

        //Viscous_Term_ZX
        viscous_components[0]= Viscous_Component_ZX(velocity_x, velocity_z,
                                                    temperature, Reynolds,
                                                    dx, dz,
                                                    i, j, k);

        //Viscous Term ZY
        viscous_components[1]= Viscous_Component_ZY(velocity_y, velocity_z,
                                                    temperature, Reynolds,
                                                    dy,dz,
                                                    i, j, k);


        //Viscous Term ZZ
        viscous_components[2]= Viscous_Component_ZZ(velocity_x, velocity_y,
                                                    velocity_z,
                                                    temperature, Reynolds,
                                                    dx, dy, dz,
                                                    i, j, k);


        viscous_total=0.;
        for(int index=0; index<3; index++)
          viscous_total += viscous_components[index];


        residual_z[k][j][i] = -convection + viscous_total + source;

      }
    }
  }


}
