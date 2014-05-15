/*  Last Modified Time-stamp: <2014-04-14 15:08:18 mike_georgiou> */
#include "../Header_Files/Data.h"
#include "../Header_Files/Macros.h"


void Flux_Evaluation_Z(double ***Flux_New_Z, double ***Speed_Int_Z,
                       double ***Rho_New, double ***Pressure_New,
                       double dz,int ldz, int ldy, int ldx)

{

  for (int  k=0; k<ldz; k++){
    for (int  j=0; j<ldy; j++){
      for (int  i=0; i<ldx; i++){

        Flux_New_Z[k][j][i] =
          9./16.*(Rho_New[k][j][i]*Speed_Int_Z[k][j][i]+
                  Rho_New[k-1][j][i]*Speed_Int_Z[k-1][j][i]) -

          1./16.*(Rho_New[k+1][j][i]*Speed_Int_Z[k+1][j][i]+
                  Rho_New[k-2][j][i]*Speed_Int_Z[k-2][j][i]) -

          dt*(9./(8.*dz)*(Pressure_New[k][j][i]-Pressure_New[k-1][j][i])-
              1./(24.*dz)*(Pressure_New[k+1][j][i]-Pressure_New[k-2][j][i]));

      }
    }
  }

}
