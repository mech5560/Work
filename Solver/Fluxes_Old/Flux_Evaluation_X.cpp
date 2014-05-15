/*  Last Modified Time-stamp: <2014-03-26 13:32:20 mike_georgiou> */
#include "../Header_Files/Data.h"
#include "../Header_Files/Macros.h"


void Flux_Evaluation_X(double ***Flux_New_X, double ***Speed_Int_X,
                       double ***Rho_New, double ***Pressure_New,
                       double dx,int ldz, int ldy, int ldx)

{

  for (int  k=0; k<ldz;   k++)
    {
      for ( int  j=0; j<ldy;  j++)
        {
          for (int   i=0; i<ldx;  i++)
            {

              Flux_New_X[k][j][i] =

                9./16.*(Rho_New[k][j][i]*Speed_Int_X[k][j][i]+
                        Rho_New[k][j][i-1]*Speed_Int_X[k][j][i-1]) -

                1./16.*(Rho_New[k][j][i+1]*Speed_Int_X[k][j][i+1]+
                        Rho_New[k][j][i-2]*Speed_Int_X[k][j][i-2]) -

                dt*(9./(8.*dx)*(Pressure_New[k][j][i]-Pressure_New[k][j][i-1])-
                    1./(24.*dx)*(Pressure_New[k][j][i+1]-Pressure_New[k][j][i-2]));

            }
        }

    }

}
