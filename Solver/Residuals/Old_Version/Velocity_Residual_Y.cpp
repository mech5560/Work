/*  Last Modified Time-stamp: <2014-04-16 12:23:50 mike_georgiou> */
/*

  This function will calculate the Velocity Residual in the Y-Momentum
  At a later stage I ll try to modify this function in order to
  reuse the calculated components of the residuals
  I am not doing that now because I am not sure that it
  will be efficient.

*/

#include"../../Header_Files/Macros.h"
#include"../../Header_Files/Data.h"




void Velocity_Residual_Y( double ***Residual_Y, double ***Speed_Z, double ***Speed_Y,
                          double ***Speed_X,double ***Delta_Y,   double ***Flux_Z,
                          double ***Flux_Y, double ***Flux_X, double ***Temperature,
                          double dz, double dx, int ldz, int ldy, int ldx)
{


  /*Variable for the calculation of the convective terms*/
  double  Res_Conv_Y[3], Convective_XY[2],Convective_ZY[2],Res_Conv_Total;


  /* Variables for the calculation of the viscous terms*/

  double Derivative_XX[4][2],Derivative_XY[4][2],
    Derivative_YY[2],Derivative_YX[2][4],Derivative_YZ[2][4],
    Derivative_ZZ[4][2], Derivative_ZY[4][2],


    Res_Visc_Y[3][4],Res_Visc_Y_Total[3], Res_Visc_Total;




  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){


        /* Convective Fluxes */


        /*X-Direction - 4th order accurate*/

        /* \delta_1 */
        Convective_XY[0]= 0.5
          *( Flux_X[k][j][i+1]*(Speed_Y[k][j][i]+Speed_Y[k][j][i+1])
             -Flux_X[k][j][i]*(Speed_Y[k][j][i]+Speed_Y[k][j][i-1]));

        /* \delta_3 */
        Convective_XY[1] = 0.5
          *(Flux_X[k][j][i+2]*(Speed_Y[k][j][i]+Speed_Y[k][j][i+3])
            -Flux_X[k][j][i-1]*(Speed_Y[k][j][i-3]+Speed_Y[k][j][i]));

        Res_Conv_Y[0] = 9./(8.*dx)*Convective_XY[0] - 1./(24.*dx)*Convective_XY[1];



        /*Y-Direction - 2nd order accurate - Non uniform grid*/
        Res_Conv_Y[1] = 0.5*
          ( (Speed_Y[k][j+1][i]+Speed_Y[k][j][i])*Flux_Y[k][j+1][i]
            -(Speed_Y[k][j-1][i]+Speed_Y[k][j][i])*Flux_Y[k][j][i])
          / (2.*Delta_Y[k][j][i]);



        /* Z-Dimension - 4th  Order accurate*/

        /* \delta_1 */
        Convective_ZY[0]= 0.5*( Flux_Z[k+1][j][i]*(Speed_Y[k][j][i]+Speed_Y[k+1][j][i])
                                -Flux_Z[k][j][i]*(Speed_Y[k][j][i]+Speed_Y[k-1][j][i]) );

        /* \delta_3 */
        Convective_ZY[1] = 0.5*(Flux_Z[k+2][j][i]*(Speed_Y[k][j][i]+Speed_Y[k+3][j][i])
                                -Flux_Z[k-1][j][i]*(Speed_Y[k][j][i]+Speed_Y[k-3][j][i]));


        Res_Conv_Y[2] = 9./(8.*dz)*Convective_ZY[0] - 1./(24.*dz)*Convective_ZY[1];



        /*Summing the Components of the Convective fluxes */
        Res_Conv_Total=0;
        for (int index=0; index<3; index++)
          Res_Conv_Total += Res_Conv_Y[index];



        /* Viscous Fluxes
           - For the calculation of these fluxes 4th order accurate

           - I derived the formulas for the fourth order accurate( with smaller stencil) scheme in the
           Intermediate-Velocity document

           - Also I need to find a proper formula for the calculation of the viscosity for liquid water*/


        /*X- Component j=1 d/dx*/

        /*Calculation of the (dv/dx) component*/

        /* Calculation of the \delta_1 component of the first derivative in the X-Direction*/
        Derivative_XX[0][0] = 9./(8.*dx)*(Speed_Y[k][j][i-1] - Speed_Y[k][j][i-2]);
        Derivative_XX[1][0] = 9./(8.*dx)*(Speed_Y[k][j][i] - Speed_Y[k][j][i-1]);
        Derivative_XX[2][0] = 9./(8.*dx)*(Speed_Y[k][j][i+1] - Speed_Y[k][j][i]);
        Derivative_XX[3][0] = 9./(8.*dx)*(Speed_Y[k][j][i+2] - Speed_Y[k][j][i+1]);


        /* Calculation of the \delta_3 component of the first derivative in the X-Direction*/
        Derivative_XX[0][1] = 1./(24.*dx)*(Speed_Y[k][j][i] - Speed_Y[k][j][i-3]);
        Derivative_XX[1][1] = 1./(24.*dx)*(Speed_Y[k][j][i+1] - Speed_Y[k][j][i-2]);
        Derivative_XX[2][1] = 1./(24.*dx)*(Speed_Y[k][j][i+2] - Speed_Y[k][j][i-1]);
        Derivative_XX[3][1] = 1./(24.*dx)*(Speed_Y[k][j][i+3] - Speed_Y[k][j][i]);



        /*Calculation of the (du/dy) component

          Each derivative will be interpolated
          in the X-direction in order to derive the necessary quantities
          For this reason I am using the 2D Derivative_XY array to  compute all the interpolated
          quantities and then combine them to get the desired result

          Each row of the Derivative_XY array will contain the two quantities needed to interpolated
          in order to obtain the desired derivative*/

        /* i-3/2*/

        /*i-3*/
        Derivative_XY[0][0] = 0.5*(Speed_X[k][j+1][i-3] - Speed_X[k][j-1][i-3])
          /(2.*Delta_Y[k][j][i-3]) ;
        /*i*/
        Derivative_XY[0][1] = 0.5*(Speed_X[k][j+1][i] - Speed_X[k][j-1][i])
          /(2.*Delta_Y[k][j][i]) ;

        /* i-1/2*/

        /*i-1*/
        Derivative_XY[1][0] = 0.5*(Speed_X[k][j+1][i-1] - Speed_X[k][j-1][i-1])
          /(2.*Delta_Y[k][j][i-1]) ;
        /*i*/
        Derivative_XY[1][1] = Derivative_XY[0][1];

        /* i+1/2*/

        /*i*/
        Derivative_XY[2][0] = Derivative_XY[0][1];
        /*i+1*/
        Derivative_XY[2][1] = 0.5*(Speed_X[k][j+1][i+1] - Speed_X[k][j-1][i+1])
          /(2.*Delta_Y[k][j][i+1]) ;

        /* i+3/2*/

        /*i*/
        Derivative_XY[3][0] = Derivative_XY[0][1];
        /*i+3*/
        Derivative_XY[3][1] = 0.5*(Speed_X[k][j+1][i+3] - Speed_X[k][j-1][i+3])
          /(2.*Delta_Y[k][j][i+3]) ;


        /* Summing the X-component of the Viscous Term of the Y-Momentum equation*/

        /* \sigma_ij i+1/2*/
        Res_Visc_Y[0][2] = /*Viscosity(i+1/2)*/

          Viscosity((Temperature[k][j][i+1]+ Temperature[k][j][i])/2. )*
          (  (Derivative_XX[2][0] - Derivative_XX[2][1])
             +  0.5*(Derivative_XY[2][1] + Derivative_XY[2][0]));


        /* \sigma_ij i-1/2*/
        Res_Visc_Y[0][1] = /*Viscosity(i-1/2)*/

          Viscosity((Temperature[k][j][i]+ Temperature[k][j][i-1])/2. )*
          ( (Derivative_XX[1][0] - Derivative_XX[1][1])
            + 0.5*( (Derivative_XY[1][1] + Derivative_XY[1][0]) ));


        /* \sigma_ij i+3/2*/
        Res_Visc_Y[0][3] = /*Viscosity(i+3/2)*/

          Viscosity((Temperature[k][j][i+1]+ Temperature[k][j][i+2])/2. )*
          (  (Derivative_XX[3][0] - Derivative_XX[3][1])
             + 0.5*( (Derivative_XY[3][1] + Derivative_XY[3][0])) );


        /* \sigma_ij i-3/2*/
        Res_Visc_Y[0][0] =/*Viscosity(i-3/2)*/

          Viscosity((Temperature[k][j][i-2]+ Temperature[k][j][i-1])/2. )*
          ( (Derivative_XX[0][0] - Derivative_XX[0][1])
            + 0.5*( (Derivative_XY[0][1] + Derivative_XY[0][0])) ) ;



        Res_Visc_Y_Total[0] =  9./(8.*dx)*(Res_Visc_Y[0][2] - Res_Visc_Y[0][1])
          -1./(24.*dx)*( Res_Visc_Y[0][3] - Res_Visc_Y[0][0]);



        /* Calculation of the Y-Component j=2  d/dy*/

        /* Calculation of the d/dy(dv/dy)*/

        /* j+1/2*/
        Derivative_YY[1] = (Speed_Y[k][j+1][i] -Speed_Y[k][j][i])
          /((Delta_Y[k][j+1][i] +Delta_Y[k][j][i]));

        /* j-1/2*/
        Derivative_YY[0] = (Speed_Y[k][j][i] -Speed_Y[k][j-1][i])
          /((Delta_Y[k][j][i] +Delta_Y[k][j-1][i]));


        /* Calculation of the d/dy(du/dx) */
        /* j-1/2*/

        /* j-1*/
        Derivative_YX[0][0] = 9./(16.*dx)*(Speed_X[k][j-1][i+1] - Speed_X[k][j-1][i-1]);
        Derivative_YX[0][1] = 1./(48.*dx)*(Speed_X[k][j-1][i+3] - Speed_X[k][j-1][i-3]);

        /* j*/
        Derivative_YX[0][2] = 9./(16.*dx)*(Speed_X[k][j][i+1] - Speed_X[k][j][i-1]);
        Derivative_YX[0][3] = 1./(48.*dx)*(Speed_X[k][j][i+3] - Speed_X[k][j][i-3]);


        /* j+1/2*/

        /* j*/
        Derivative_YX[1][0] = Derivative_YX[0][2];
        Derivative_YX[1][1] = Derivative_YX[0][3];

        /* j+1*/
        Derivative_YX[1][2] = 9./(16.*dx)*(Speed_X[k][j+1][i+1] - Speed_X[k][j+1][i-1]);
        Derivative_YX[1][3] = 1./(48.*dx)*(Speed_X[k][j+1][i+3] - Speed_X[k][j+1][i-3]);



        /* Calculation of the d/dy(dw/dz), i am using the \delta_1 -\delta_3 scheme
           maybe i could use the \delta_1 -\delta_2 scheme*/


        /* j-1/2*/

        /* j-1*/
        Derivative_YZ[0][0] = 9./(16.*dz)*(Speed_Z[k+1][j-1][i] - Speed_Z[k-1][j-1][i]);
        Derivative_YZ[0][1] = 1./(48.*dz)*(Speed_Z[k+3][j-1][i] - Speed_Z[k-3][j-1][i]);

        /* j*/
        Derivative_YZ[0][2] = 9./(16.*dz)*(Speed_Z[k+1][j][i] - Speed_Z[k-1][j][i]);
        Derivative_YZ[0][3] = 1./(48.*dz)*(Speed_Z[k+3][j][i] - Speed_Z[k-3][j][i]);


        /* j+1/2*/

        /* j*/
        Derivative_YZ[1][0] = Derivative_YZ[0][2];
        Derivative_YZ[1][1] = Derivative_YZ[0][3];

        /* j+1*/
        Derivative_YZ[1][2] = 9./(16.*dz)*(Speed_Z[k+1][j+1][i] - Speed_Z[k-1][j+1][i]);
        Derivative_YZ[1][3] = 1./(48.*dz)*(Speed_Z[k+3][j+1][i] - Speed_Z[k-3][j+1][i]);



        /* Calculating the components of the total residular Res_Visc_Total[1]*/

        /* \sigma_ij j+1/2*/

        Res_Visc_Y[1][1] = /*Viscosity_j+1/2*/

          Viscosity((Temperature[k][j][i]+ Temperature[k][j][i])/2. )*
          ( 4./3.*Derivative_YY[1]
            -2./6.*( (Derivative_YX[1][2] - Derivative_YX[1][3])
                     +   (Derivative_YX[1][0] - Derivative_YX[1][1])

                     +   (Derivative_YZ[1][2] - Derivative_YZ[1][3])
                     +   (Derivative_YZ[1][0] - Derivative_YZ[1][1]) ) );


        /* \sigma_ij j-1/2*/

        Res_Visc_Y[1][0] = /* \viscosity j-1/2*/

          Viscosity((Temperature[k][j][i]+ Temperature[k][j][i])/2. )*
          (4./3.*Derivative_YY[0]
           -2./6.*( (Derivative_YX[0][2] - Derivative_YX[0][3])
                    +  (Derivative_YX[0][0] - Derivative_YX[0][1])
                    +  (Derivative_YZ[0][2] - Derivative_YZ[0][3])
                    +  (Derivative_YZ[0][0] - Derivative_YZ[0][1])) );


        Res_Visc_Y_Total[1] = 1./(Reynolds * (2.* Delta_Y[k][j][i]))
          *(Res_Visc_Y[1][1] - Res_Visc_Y[1][0]);



        /* Calculation of the Z-Component j=3*/

        /* Calculation of the d/dz(dv/dz) component*/


        /* Calculation of the \delta_1 component of the first derivative in the Z-Direction*/
        Derivative_ZZ[0][0] = 9./(8.*dz)*(Speed_Y[k-1][j][i] - Speed_Y[k-2][j][i]);
        Derivative_ZZ[1][0] = 9./(8.*dz)*(Speed_Y[k][j][i] - Speed_Y[k-1][j][i]);
        Derivative_ZZ[2][0] = 9./(8.*dz)*(Speed_Y[k+1][j][i] - Speed_Y[k][j][i]);
        Derivative_ZZ[3][0] = 9./(8.*dz)*(Speed_Y[k+2][j][i] - Speed_Y[k+1][j][i]);


        /* Calculation of the \delta_3 component of the first derivative in the Z-Direction*/
        Derivative_ZZ[0][1] = 1./(24.*dz)*(Speed_Y[k][j][i] - Speed_Y[k-3][j][i]);
        Derivative_ZZ[1][1] = 1./(24.*dz)*(Speed_Y[k+1][j][i] - Speed_Y[k-2][j][i]);
        Derivative_ZZ[2][1] = 1./(24.*dz)*(Speed_Y[k+2][j][i] - Speed_Y[k-1][j][i]);
        Derivative_ZZ[3][1] = 1./(24.*dz)*(Speed_Y[k+3][j][i] - Speed_Y[k][j][i]);




        /* Calculation of the d/dz(dw/dy) component*/

        /* k-3/2*/

        /*k-3*/
        Derivative_ZY[0][0] = 0.5*(Speed_Z[k-3][j+1][i] - Speed_Z[k-3][j-1][i])
          /(2.*Delta_Y[k-3][j][i]);

        /*k*/
        Derivative_ZY[0][1] = 0.5*(Speed_Z[k][j+1][i] - Speed_Z[k][j-1][i])
          /(2.*Delta_Y[k][j][i]);

        /* k-1/2*/

        /*k-1*/
        Derivative_ZY[1][0] = 0.5*(Speed_Z[k-1][j+1][i] - Speed_Z[k-1][j-1][i])
          /(2.*Delta_Y[k-1][j][i]);

        /*k*/
        Derivative_ZY[1][1] = Derivative_ZY[0][1];


        /* k+1/2*/

        /*k*/
        Derivative_ZY[2][0] = Derivative_ZY[0][1];
        /*k+1*/
        Derivative_ZY[2][1] = 0.5*(Speed_Z[k+1][j+1][i] - Speed_Z[k+1][j-1][i])
          /(2.*Delta_Y[k+1][j][i]);

        /* k+3/2*/
        /*k*/
        Derivative_ZY[3][0] = Derivative_ZY[0][1];
        /*k+3*/
        Derivative_ZY[3][1] = 0.5*(Speed_Z[k+3][j+1][i] - Speed_Z[k+3][j-1][i])
          /(2.*Delta_Y[k+3][j][i]);



        /* Calculating the components of the total residular Res_Visc_Y_Total[2]*/



        /* \sigma_ij k+1/2*/

        Res_Visc_Y[2][2] = /*Viscosity_k+1/2*/

          Viscosity((Temperature[k+1][j][i]+ Temperature[k][j][i])/2. )*
          ((Derivative_ZZ[2][0]- Derivative_ZZ[2][1])
           +0.5*( Derivative_ZY[2][1] + Derivative_ZY[2][0] ));

        /* \sigma_ij k-1/2*/

        Res_Visc_Y[2][1] = /*Viscosity_k-1/2*/

          Viscosity((Temperature[k][j][i]+ Temperature[k-1][j][i])/2 )*
          ( (Derivative_ZZ[1][0]- Derivative_ZZ[1][1])
            +0.5*( Derivative_ZY[1][1] + Derivative_ZY[1][0] ));


        /* \sigma_ij k+3/2*/
        Res_Visc_Y[2][3] = /*Viscosity_k+3/2*/

          Viscosity((Temperature[k+2][j][i]+ Temperature[k+1][j][i])/2 )*
          ((Derivative_ZZ[3][0]- Derivative_ZZ[3][1])
           +0.5*( Derivative_ZY[3][1] + Derivative_ZY[3][0] ));


        /* \sigma_ij k-3/2*/
        Res_Visc_Y[2][0] = /*Viscosity_k-3/2*/
          Viscosity((Temperature[k-1][j][i]+ Temperature[k-2][j][i])/2 )*
          ((Derivative_ZZ[0][0]- Derivative_ZZ[0][1])
           +0.5*( Derivative_ZY[0][1] + Derivative_ZY[0][0] ));



        Res_Visc_Y_Total[2] = 1/Reynolds
          *( 9./(8.*dz)*(Res_Visc_Y[2][2] -Res_Visc_Y[2][1])
             -1./(24.*dz)*(Res_Visc_Y[2][3] -Res_Visc_Y[2][0]));



        Res_Visc_Total=0.;
        for(int index=0; index<3; index++)
          Res_Visc_Total+=Res_Visc_Y_Total[index];


        Residual_Y[k][j][i] = -Res_Conv_Total + Res_Visc_Total;

      }
    }
  }


}
