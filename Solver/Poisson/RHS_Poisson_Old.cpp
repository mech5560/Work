/*  Last Modified Time-stamp: <2014-04-01 17:41:08 mike_georgiou> */ 
/* Right hand side of the poisson equation

   - In this function, the B vector of the $Ax = B$ equation will be defined.

   *********** Notes ****************

   - Dont forget that triple pointers (***) are faster than
   the mapped 1d array.

   -I already have the relations for the coefficients of the RHS of the Poisson eqn.
   -They are the same with the coefficients of the diffussion term but divided with dx and not dx^2

*/


#include "./../../Header_Files/Data.h"
#include "./../../Header_Files/Macros.h"


 void RHS_Poisson_Old( double *RHS, double ***Speed_Int_X,
		       double ***Speed_Int_Y, double ***Speed_Int_Z,
		       double ***Rho_Star, double ***Rho_N, double ***Rho_Nm1, double ***Delta_Y,
		       double dz, double dx,int ldz, int ldy, int ldx)
{


  double Density_Der,Auxiliary_Der[3];

  double Coefficients_Z[2], Coefficients_X[2];


  Coefficients_X[0]= 9./(16.*dx);
  Coefficients_X[1]= 1./(48.*dx);

  Coefficients_Z[0]= 9./(16.*dz);
  Coefficients_Z[1]= 1./(48.*dz);



  for (int k=0; k<ldz; k++)
    {
      for (int j=0; j<ldy; j++)
        {
          for (int i=0; i<ldx; i++)
            {



              /* Calculation of the time derivative of the Density */
              Density_Der = 0.5*( 3.*Rho_Star[k][j][i] - 4.*Rho_N[k][j][i] + Rho_Nm1[k][j][i] )
                /(dt*dt);


              /* Calculation of the Auxiliary Derivatives in the X-Direction*/

              Auxiliary_Der[0] =
                (Coefficients_X[0]*(Rho_Star[k][j][i+1]*Speed_Int_X[k][j][i+1] -
				    Rho_Star[k][j][i-1]* Speed_Int_X[k][j][i-1]) -
                
		 Coefficients_X[1]*(Rho_Star[k][j][i+3]*Speed_Int_X[k][j][i+3] -
				    Rho_Star[k][j][i-3]*Speed_Int_X[k][j][i-3] ))/dt;


              /* Calculation of the Auxiliary Derivatives in the Y-Direction*/

	      Auxiliary_Der[1] = (0.5*(Rho_Star[k][j+1][i]*Speed_Int_Y[k][j+1][i]-
				       Rho_Star[k][j-1][i]*Speed_Int_Y[k][j-1][i])
				  /(2.*Delta_Y[k][j][i]))/dt;


              /* Calculation of the Auxiliary Derivatives in the Z-Direction*/

              Auxiliary_Der[2] =
                (Coefficients_Z[0]*(Rho_Star[k+1][j][i]*Speed_Int_Z[k+1][j][i] -
                                   Rho_Star[k-1][j][i]*Speed_Int_Z[k-1][j][i]) -
                Coefficients_Z[1]*(Rho_Star[k+3][j][i]*Speed_Int_Z[k+3][j][i] -
                                   Rho_Star[k-3][j][i]*Speed_Int_Z[k-3][j][i] ))/dt;

              /* A(k,j,i) Maps the location of the current
		 position to the one-dimensional Vector */

              RHS[A(k,j,i)+1] =
                Auxiliary_Der[0]+
                Auxiliary_Der[1]+
                Auxiliary_Der[2]+
                Density_Der;

            }
        }
    }

}
