/*  Last Modified Time-stamp: <2014-04-15 14:14:13 mike_georgiou> */
#include "../Header_Files/Data.h"
#include "../Header_Files/Macros.h"




void Flux_Evaluation_Y(double ***Flux_New_Y, double ***Speed_Int_Y,
                       double ***Rho_New, double ***Pressure_New,
                       double ***Points_Y,int ldz, int ldy, int ldx)

{

  for (int k=0; k<ldz;  k++){
    for (int j=1; j<ldy-1;  j++){
      for (int i=0; i<ldx; i++){

        Flux_New_Y[k][j][i] =
					
          0.5*(Rho_New[k][j][i]*Speed_Int_Y[k][j][i]+
               Rho_New[k][j-1][i]*Speed_Int_Y[k][j-1][i]) -

          dt*(Pressure_New[k][j][i]-Pressure_New[k][j-1][i])/
          (Points_Y[k][j][i]+Points_Y[k][j-1][i]);


        /*Forcing the Wall BC for my case */
        Flux_New_Y[k][0][i]=0.;
        Flux_New_Y[k][ldy-1][i]=0.;

      }
    }
  }

}
