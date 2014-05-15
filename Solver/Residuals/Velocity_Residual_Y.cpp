/*  Last Modified Time-stamp: <2014-05-15 12:48:02 mike_georgiou> */
/*

  This function will calculate the Velocity Residual in the Y-Momentum
  At a later stage I ll try to modify this function in order to
  reuse the calculated components of the residuals
  I am not doing that now because I am not sure that it
  will be efficient.

*/

#include"Residuals.h"
#include"Residuals-inl.h"

#include"../../Header_Files/Data.h"

void Velocity_Residual_Y( double*** residual_y,double*** velocity_x,
                          double*** velocity_y,double*** velocity_z,
                          double*** flux_x,double*** flux_y, double***flux_z,
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


        //Convection Term
        convection= Convection_Term(velocity_y,
                                    flux_x, flux_y,flux_z,
                                    dx, dy, dz,
                                    i, j, k);

        //Viscous Term YX
        viscous_components[0]= Viscous_Component_YX(velocity_x, velocity_y,
                                                    temperature, Reynolds,
                                                    dx, dy,
                                                    i, j, k);

        //Viscous Term YY
        viscous_components[1]= Viscous_Component_YY(velocity_x, velocity_y,
                                                    velocity_z,
                                                    temperature, Reynolds,
                                                    dx, dy,dz,
                                                    i, j, k);


        //Viscous Term YZ
        viscous_components[2]= Viscous_Component_YZ(velocity_y, velocity_z,
                                                    temperature, Reynolds,
                                                    dy,dz,
                                                    i, j, k);

        viscous_total=0.;
        for(int index=0; index<3; index++)
          viscous_total += viscous_components[index];

        residual_y[k][j][i] = -convection + viscous_total + source;


      }
    }
  }


}
