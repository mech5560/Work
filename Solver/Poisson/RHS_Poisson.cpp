/*Last Modified Time-stamp: <2014-04-25 13:15:06 mike_georgiou>

  Author: Michail Georgiou. Date Created: Tue Apr  1 17:00:03 CEST 2014

  Description:

  In this version of the RHS_Poisson function I am using the
  scheme that is presented in the Lessani-Papalexandris paper eq.66

  One potential difference between this scheme and the one that i am using
  is that this scheme is conserving the quantities within each cell in comparison
  with the previous one whis is expanded to greater number of cells.

  It would be nice to compare the effectiveness of these two versions of the code
*/

#include "./../../Header_Files/Data.h"
#include "./../../Header_Files/Macros.h"


void RHS_Poisson( double *RHS, double ***Speed_Int_X,
                  double ***Speed_Int_Y, double ***Speed_Int_Z,
                  double ***Rho_Star, double ***Rho_N, double ***Rho_Nm1, double ***Delta_Y,
                  double dz, double dx,int ldz, int ldy, int ldx)
{

  double Density_Der,Auxiliary_Der[3], Derivative_X[4][2], Derivative_Z[4][2];

  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

        //Calculation of the time derivative of the Density

        // Density_Der = (3.0*Rho_Star[k][j][i]-4.0*Rho_N[k][j][i]+
        // 1.0*Rho_Nm1[k][j][i])/(2.*dt);

        //First Order, Just for Checking
        Density_Der = (Rho_Star[k][j][i]-Rho_N[k][j][i])/(dt);


        /*calculation of the Auxiliary Derivatives in the X-Direction*/

        /*\delta_1 int_1 i-1/2*/
        Derivative_X[0][0] = (9./16.*(Rho_Star[k][j][i]*Speed_Int_X[k][j][i] +
                                      Rho_Star[k][j][i-1]*Speed_Int_X[k][j][i-1]));

        /*\delta_1 int_1 i+1/2*/
        Derivative_X[0][1] = (9./16.*(Rho_Star[k][j][i+1]*Speed_Int_X[k][j][i+1] +
                                      Rho_Star[k][j][i]*Speed_Int_X[k][j][i]));

        /*\delta_1 int_3 i-1/2*/
        Derivative_X[1][0] = (1./16.*(Rho_Star[k][j][i+1]*Speed_Int_X[k][j][i+1] +
                                      Rho_Star[k][j][i-2]*Speed_Int_X[k][j][i-2]));

        /*\delta_1 int_3 i+1/2*/
        Derivative_X[1][1] = (1./16.*(Rho_Star[k][j][i+2]*Speed_Int_X[k][j][i+2] +
                                      Rho_Star[k][j][i-1]*Speed_Int_X[k][j][i-1]));


        /*\delta_3 int_1 i-3/2*/
        Derivative_X[2][0] = (9./16.*(Rho_Star[k][j][i-1]*Speed_Int_X[k][j][i-1] +
                                      Rho_Star[k][j][i-2]*Speed_Int_X[k][j][i-2]));

        /*\delta_3 int_1 i+3/2*/
        Derivative_X[2][1] = (9./16.*(Rho_Star[k][j][i+2]*Speed_Int_X[k][j][i+2] +
                                      Rho_Star[k][j][i+1]*Speed_Int_X[k][j][i+1]));


        /*\delta_3 int_3 i-3/2*/
        Derivative_X[3][0] = (1./16.*(Rho_Star[k][j][i]*Speed_Int_X[k][j][i] +
                                      Rho_Star[k][j][i-3]*Speed_Int_X[k][j][i-3]));

        /*\delta_3 int_3 i+3/2*/
        Derivative_X[3][1] = (1./16.*(Rho_Star[k][j][i+3]*Speed_Int_X[k][j][i+3] +
                                      Rho_Star[k][j][i]*Speed_Int_X[k][j][i]));


        Auxiliary_Der[0] =
          9./(8.*dx)*((Derivative_X[0][1]-Derivative_X[0][0])
                      -(Derivative_X[1][1]-Derivative_X[1][0]))
          -1./(24.*dx)*((Derivative_X[2][1]-Derivative_X[2][0])
                        -(Derivative_X[3][1]-Derivative_X[3][0]));



        /* Calculation of the Auxiliary Derivatives in the Y-Direction*/

        Auxiliary_Der[1] = (0.5*(Rho_Star[k][j+1][i]*Speed_Int_Y[k][j+1][i]-
                                 Rho_Star[k][j-1][i]*Speed_Int_Y[k][j-1][i])
                            /(2.*Delta_Y[k][j][i]));



        /* Calculation of the Auxiliary Derivatives in the Z-Direction*/

        /*\delta_1 int_1 i-1/2*/
        Derivative_Z[0][0] = (9./16.*(Rho_Star[k][j][i]*Speed_Int_Z[k][j][i] +
                                      Rho_Star[k-1][j][i]*Speed_Int_Z[k-1][j][i]));

        /*\delta_1 int_1 i+1/2*/
        Derivative_Z[0][1] = (9./16.*(Rho_Star[k+1][j][i]*Speed_Int_Z[k+1][j][i] +
                                      Rho_Star[k][j][i]*Speed_Int_Z[k][j][i]));


        /*\delta_1 int_3 i-1/2*/
        Derivative_Z[1][0] = (1./16.*(Rho_Star[k+1][j][i]*Speed_Int_Z[k+1][j][i]+
                                      Rho_Star[k-2][j][i]*Speed_Int_Z[k-2][j][i]));

        /*\delta_1 int_3 i+1/2*/
        Derivative_Z[1][1] = (1./16.*(Rho_Star[k+2][j][i]*Speed_Int_Z[k+2][j][i] +
                                      Rho_Star[k-1][j][i]*Speed_Int_Z[k-1][j][i]));


        /*\delta_3 int_1 i-3/2*/
        Derivative_Z[2][0] = (9./16.*(Rho_Star[k-1][j][i]*Speed_Int_Z[k-1][j][i] +
                                      Rho_Star[k-2][j][i]*Speed_Int_Z[k-2][j][i]));

        /*\delta_3 int_1 i+3/2*/
        Derivative_Z[2][1] = (9./16.*(Rho_Star[k+2][j][i]*Speed_Int_Z[k+2][j][i] +
                                      Rho_Star[k+1][j][i]*Speed_Int_Z[k+1][j][i]));



        /*\delta_3 int_3 i-3/2*/
        Derivative_Z[3][0] = (1./16.*(Rho_Star[k][j][i]*Speed_Int_Z[k][j][i] +
                                      Rho_Star[k-3][j][i]*Speed_Int_Z[k-3][j][i]));

        /*\delta_3 int_3 i+3/2*/
        Derivative_Z[3][1] = (1./16.*(Rho_Star[k+3][j][i]*Speed_Int_Z[k+3][j][i] +
                                      Rho_Star[k][j][i]*Speed_Int_Z[k][j][i]));


        /*Summing the Auziliary_Derivative components*/
        Auxiliary_Der[2] =
          9./(8.*dz)*((Derivative_Z[0][1]-Derivative_Z[0][0])
                      -(Derivative_Z[1][1]-Derivative_Z[1][0]))
          -1./(24.*dz)*((Derivative_Z[2][1]-Derivative_Z[2][0])
                        -(Derivative_Z[3][1]-Derivative_Z[3][0]));

        double Deriv_Sum=0.;
        for (int i=0; i<3; i++)
          Deriv_Sum += Auxiliary_Der[i];



        /* Calculation of the Auxiliary Derivatives in the Z-Direction*/


        /* A(k,j,i) Maps the location of the current
           position to the one-dimensional Vector */

        RHS[A(k,j,i)+1] = (Density_Der+Deriv_Sum)/dt;

      }
    }
  }

}
