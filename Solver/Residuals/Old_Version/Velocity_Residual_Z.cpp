/*  Last Modified Time-stamp: <2014-03-26 17:30:28 mike_georgiou> */ 
/*   Last Modified Time-stamp: <2014-03-19 18:41:25 mike_georgiou>  
This function will calculate the Velocity Residual in the X-Direction. 


At a later stage I ll try to modify this function in order to 
reuse the calculated components of the residuals

I am not doing that now because I am not sure that it 
will be efficient.

*/

#include"../../Header_Files/Macros.h"
#include"../../Header_Files/Data.h"


void Velocity_Residual_Z( double ***Residual_Z, double ***Speed_Z, double ***Speed_Y, 
			  double ***Speed_X,double ***Delta_Y,   double ***Flux_Z,
			  double ***Flux_Y, double ***Flux_X,double ***Temperature,
			  double dz, double dx, int ldz, int ldy, int ldx){


  /*Variable for the calculation of the convective terms*/
  double  Res_Conv_Z[3], Convective_ZZ[2],Convective_ZX[2],Res_Conv_Total;
   

  /* Variables for the calculation of the viscous terms*/
  
  double Derivative_XX[4][2],Derivative_XZ[4][4],
         Derivative_YY[2],Derivative_YZ[2][4],
         Derivative_ZZ[4][2], Derivative_ZX[4][4],Derivative_ZY[4][2],
         Res_Visc_Z[3][4],Res_Visc_Z_Total[3], Res_Visc_Total;


  

  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

	

	/* Convective Fluxes */


        /*X-Direction - 4th order accurate*/

	/* \delta_1 */
	Convective_ZX[0]= 0.5*(Flux_X[k][j][i+1]*
			       (Speed_Z[k][j][i+1]+Speed_Z[k][j][i]) 
			     - Flux_X[k][j][i]*
			       (Speed_Z[k][j][i-1]+Speed_Z[k][j][i]) );

	/* \delta_3 */
	Convective_ZX[1] = 0.5*(Flux_X[k][j][i+2]*
				(Speed_Z[k][j][i+3]+Speed_Z[k][j][i])
				-Flux_X[k][j][i-1]*
				(Speed_Z[k][j][i-3]+Speed_Z[k][j][i]) );

	Res_Conv_Z[0] = 9./(8.*dx)*Convective_ZX[0] - 1./(24.*dx)*Convective_ZX[1]; 
		      


	/*Y-Direction - 2nd order accurate - Non uniform grid*/

	    Res_Conv_Z[1] = 
	      0.5*( (Speed_Z[k][j+1][i]+Speed_Z[k][j][i])*Flux_Y[k][j+1][i] 
		    -(Speed_Z[k][j-1][i]+Speed_Z[k][j][i])*Flux_Y[k][j][i]) 
	      / (2.*Delta_Y[k][j][i] );


        /* Z-Dimension - 4th  Order accurate*/

	/* \delta_1 */
	Convective_ZZ[0]= 0.5*(Flux_Z[k+1][j][i]*
			       (Speed_Z[k+1][j][i]+Speed_Z[k][j][i]) 
			       - Flux_Z[k][j][i]*
			       (Speed_Z[k-1][j][i]+Speed_Z[k][j][i]));

	/* \delta_3 */
	Convective_ZZ[1] = 0.5*(Flux_Z[k+2][j][i]*
				(Speed_Z[k+3][j][i]+Speed_Z[k][j][i])
				-Flux_Z[k-1][j][i]*
				(Speed_Z[k-3][j][i]+Speed_Z[k][j][i]) );


	Res_Conv_Z[2] = 9./(8.*dz)*Convective_ZZ[0] - 1./(24.*dz)*Convective_ZZ[1]; 



	/*Summing the Components of the Convective fluxes */
	Res_Conv_Total=0.;
	for (int index=0; index<3; index++)
	  Res_Conv_Total +=Res_Conv_Z[index];



/* Viscous Fluxes 
   
  - For the calculation of these fluxes,(4th order accurate)

  - I derived the formulas for the fourth order accurate( with smaller stencil) scheme in the 
    Intermediate-Velocity document

*/ 

	/*X- Component j=1 d/dx*/


	/* Calculation of the d/dx(dw/dx) component*/


	/* Calculation of the \delta_1 component of the first derivative in the Z-Direction*/
	Derivative_XX[0][0] = 9./(8.*dx)*(Speed_Z[k][j][i-1] - Speed_Z[k][j][i-2]);
	Derivative_XX[1][0] = 9./(8.*dx)*(Speed_Z[k][j][i] - Speed_Z[k][j][i-1]);
	Derivative_XX[2][0] = 9./(8.*dx)*(Speed_Z[k][j][i+1] - Speed_Z[k][j][i]);
	Derivative_XX[3][0] = 9./(8.*dx)*(Speed_Z[k][j][i+2] - Speed_Z[k][j][i+1]);


	/* Calculation of the \delta_3 component of the first derivative in the Z-Direction*/
	Derivative_XX[0][1] = 1./(24.*dx)*(Speed_Z[k][j][i] - Speed_Z[k][j][i-3]);
	Derivative_XX[1][1] = 1./(24.*dx)*(Speed_Z[k][j][i+1] - Speed_Z[k][j][i-2]);
	Derivative_XX[2][1] = 1./(24.*dx)*(Speed_Z[k][j][i+2] - Speed_Z[k][j][i-1]);
	Derivative_XX[3][1] = 1./(24.*dx)*(Speed_Z[k][j][i+3] - Speed_Z[k][j][i]);




	/* Calculation of the d/dx(du/dz) component*/

	/**** i-3/2 ****/

	/* i-3 */
	/*\delta_1*/
	Derivative_XZ[0][0] = 2./(3.*dz)
	  *(Speed_X[k+1][j][i-3] - Speed_X[k-1][j][i-3] );
	/*\delta_2*/
	Derivative_XZ[0][1] = 1./(12.*dz)
	  *(Speed_X[k+2][j][i-3] - Speed_X[k-2][j][i-3]);

	/* i */
	/*\delta_1*/
	Derivative_XZ[0][2] = 2./(3.*dz)
	  *(Speed_X[k+1][j][i] - Speed_X[k-1][j][i]);
	/*\delta_2*/
	Derivative_XZ[0][3] = 1./(12.*dz)
	  *(Speed_X[k+2][j][i] - Speed_X[k-2][j][i]);


	/**** i-1/2 ****/
	
	/*i-1*/
	/*\delta_1*/
	Derivative_XZ[1][0] = 2./(3.*dz)
	  *(Speed_X[k+1][j][i-1] - Speed_X[k-1][j][i-1] );
	/*\delta_2*/
	Derivative_XZ[1][1] = 1./(12.*dz)
	  *(Speed_X[k+2][j][i-1] - Speed_X[k-2][j][i-1]);

	/*i*/
	/*\delta_1*/	
	Derivative_XZ[1][2]=Derivative_XZ[0][2];
	/*\delta_2*/
	Derivative_XZ[1][3]=Derivative_XZ[0][3];

	/**** i+1/2 ****/
	
	/*i*/
	/*\delta_1*/
	Derivative_XZ[2][0] = Derivative_XZ[0][2];
	/*\delta_2*/
	Derivative_XZ[2][1] = Derivative_XZ[0][3];

	/*i+1*/
	/*\delta_1*/
	Derivative_XZ[2][2] = 2./(3.*dz)
	  *(Speed_X[k+1][j][i+1] - Speed_X[k-1][j][i+1] );
	/*\delta_2*/
	Derivative_XZ[2][3] = 1./(12.*dz)
	  *(Speed_X[k+2][j][i+1] - Speed_X[k-2][j][i+1]);



	/**** i+3/2 ****/
	
	/*i*/
	/*\delta_1*/
	Derivative_XZ[3][0] = Derivative_XZ[0][2];
	/*\delta_2*/
	Derivative_XZ[3][1] = Derivative_XZ[0][3];


	/*i+3*/
	/*\delta_1*/
	Derivative_XZ[3][2] = 2./(3.*dz)
	  *(Speed_X[k+1][j][i+3] - Speed_X[k-1][j][i+3] );
	/*\delta_2*/
	Derivative_XZ[3][3] = 1./(12.*dz)
	  *(Speed_X[k+2][j][i+3] - Speed_X[k-2][j][i+3]);




      /* Calculating the components of the total residular Res_Visc_Z_Total[0]*/

	    

	    /* \sigma_ij i+1/2*/

	Res_Visc_Z[0][2] = /*Viscosity_i+1/2*/ 
	  
	  Viscosity((Temperature[k][j][i+1]+ Temperature[k][j][i])/2. )*
	  ( Derivative_XX[2][0]- Derivative_XX[2][1]
	    +0.5*( (Derivative_XZ[2][2]-Derivative_XZ[2][3]) 
		   +(Derivative_XZ[2][0]-Derivative_XZ[2][1]) ));


	    /* \sigma_ij i-1/2*/

	    Res_Visc_Z[0][1] = /*Viscosity_i-1/2*/ 

	      Viscosity((Temperature[k][j][i-1]+ Temperature[k][j][i])/2. )*
	      ( Derivative_XX[1][0]- Derivative_XX[1][1]
		+0.5*( (Derivative_XZ[1][0]-Derivative_XZ[1][1]) 
		       +(Derivative_XZ[1][2]-Derivative_XZ[1][3]) ));


	    /* \sigma_ij i+3/2*/

	    Res_Visc_Z[0][3] = /*Viscosity_i+3/2*/ 
	       
	      Viscosity((Temperature[k][j][i+1]+ Temperature[k][j][i+2])/2. )*
	      (Derivative_XX[3][1]- Derivative_XX[3][2]
	       +0.5*( (Derivative_XZ[3][0]-Derivative_XZ[3][1]) 
		      +(Derivative_XZ[3][2]-Derivative_XZ[3][3]) ));

	    /* \sigma_ij i-3/2*/

	    Res_Visc_Z[0][0] = /*Viscosity_i-3/2*/ 

	      Viscosity((Temperature[k][j][i-1]+ Temperature[k][j][i-2])/2. )*
	      ( Derivative_XX[0][1]- Derivative_XX[0][2]
		+0.5*( (Derivative_XZ[0][0]-Derivative_XZ[0][1]) 
		       +(Derivative_XZ[0][2]-Derivative_XZ[0][3]) ));



	    Res_Visc_Z_Total[0] = 1./Reynolds
	      *( 9./(8.*dx)*(Res_Visc_Z[0][2] -Res_Visc_Z[0][1])  
		 -1./(24.*dx)*(Res_Visc_Z[0][3] -Res_Visc_Z[0][0]));



	/* Calculation of the Y-Component j=2 d/dy*/



	/* Calculation of the d/dy(dw/dy)*/

	    /* j+1/2*/
	    Derivative_YY[1] = (Speed_Z[k][j+1][i] -Speed_Z[k][j][i])
	      /((Delta_Y[k][j+1][i] +Delta_Y[k][j][i]));

	    /* j-1/2*/
	    Derivative_YY[0] = (Speed_Z[k][j][i] -Speed_Z[k][j-1][i])
	      /((Delta_Y[k][j][i] +Delta_Y[k][j-1][i]));




	/* Calculation of the d/dy(dv/dz) */


	/* j-1/2*/

	    /* j-1*/
	    Derivative_YZ[0][0] = 9./(16.*dz)*(Speed_Y[k+1][j-1][i] - Speed_Y[k-1][j-1][i]);
	    Derivative_YZ[0][1] = 1./(48.*dz)*(Speed_Y[k+3][j-1][i] - Speed_Y[k-3][j-1][i]);

	    /* j*/
	    Derivative_YZ[0][2] = 9./(16.*dz)*(Speed_Y[k+1][j][i] - Speed_Y[k-1][j][i]);
	    Derivative_YZ[0][3] = 1./(48.*dz)*(Speed_Y[k+3][j][i] - Speed_Y[k-3][j][i]);


	/* j+1/2*/

	    /* j*/
	    Derivative_YZ[1][0] = Derivative_YZ[0][2];
	    Derivative_YZ[1][1] = Derivative_YZ[0][3];

	    /* j+1*/
	    Derivative_YZ[1][2] = 9./(16.*dz)*(Speed_Y[k+1][j+1][i] - Speed_Y[k-1][j+1][i]);
	    Derivative_YZ[1][3] = 1./(48.*dz)*(Speed_Y[k+3][j+1][i] - Speed_Y[k-3][j+1][i]);




      /* Calculating the components of the total residular Res_Visc_Total[1]*/
	 
	    /* \sigma_ij j+1/2*/

	 Res_Visc_Z[1][1] = /*Viscosity_j+1/2*/ 

	 Viscosity((Temperature[k][j+1][i]+ Temperature[k][j][i])/2. )*
	   (Derivative_YY[1] + 0.5*( (Derivative_YZ[1][2] - Derivative_YZ[1][3]) 
				    +(Derivative_YZ[1][0] - Derivative_YZ[1][1]) ));
	    

	    /* \sigma_ij j-1/2*/
	    Res_Visc_Z[1][0] = /* \viscosity j-1/2*/ 

	      Viscosity((Temperature[k][j][i]+ Temperature[k][j][i])/2. )*
	       ( Derivative_YY[0]  
		 +0.5*( (Derivative_YZ[0][2] - Derivative_YZ[0][3]) 
			+(Derivative_YZ[0][0] - Derivative_YZ[0][1]) ));


	    Res_Visc_Z_Total[1] = 1./(Reynolds
				     *(2.*Delta_Y[k][j][i]))
	                             *(Res_Visc_Z[1][1] - Res_Visc_Z[1][0]);



	/* Calculation of the Z-Component j=3 d/dz*/


  	/*Calculation of the d/dz(dw/dz) component*/

	/* Calculation of the \delta_1 component of the first derivative in the Z-Direction*/
	Derivative_ZZ[0][0] = 9./(8.*dz)*(Speed_Z[k-1][j][i] - Speed_Z[k-2][j][i]);
	Derivative_ZZ[1][0] = 9./(8.*dz)*(Speed_Z[k][j][i] - Speed_Z[k-1][j][i]);
	Derivative_ZZ[2][0] = 9./(8.*dz)*(Speed_Z[k+1][j][i] - Speed_Z[k][j][i]);
	Derivative_ZZ[3][0] = 9./(8.*dz)*(Speed_Z[k+2][j][i] - Speed_Z[k+1][j][i]);


	/* Calculation of the \delta_3 component of the first derivative in the Z-Direction*/
	Derivative_ZZ[0][1] = 1./(24.*dz)*(Speed_Z[k][j][i] - Speed_Z[k-3][j][i]);
	Derivative_ZZ[1][1] = 1./(24.*dz)*(Speed_Z[k+1][j][i] - Speed_Z[k-2][j][i]);
	Derivative_ZZ[2][1] = 1./(24.*dz)*(Speed_Z[k+2][j][i] - Speed_Z[k-1][j][i]);
	Derivative_ZZ[3][1] = 1./(24.*dz)*(Speed_Z[k+3][j][i] - Speed_Z[k][j][i]);



	/*Calculation of the d/dz(dv/dy) component 
	  
	  Each derivative will be interpolated 
	  in the Z-direction in order to derive the necessary quantities
	  For this reason I am using the 2D  array, Derivative_ZY, to  compute all the interpolated
	  quantities and then combine them to get the desired result
	 
	  Each row of the Derivative_ZY array will contain the two quantities needed to interpolated
	  in order to obtain the desired derivative*/


	/* k-3/2*/
	Derivative_ZY[0][0] = 0.5*(Speed_Y[k-3][j+1][i] - Speed_Y[k-3][j-1][i])
	                      /(2.* Delta_Y[k-3][j][i] ) ;

	Derivative_ZY[0][1] = 0.5* (Speed_Y[k][j+1][i] - Speed_Y[k][j-1][i])
	                      /(2.* Delta_Y[k][j][i] ) ;

	/* k-1/2*/
	Derivative_ZY[1][0] =  0.5*(Speed_Y[k-1][j+1][i] - Speed_Y[k-1][j-1][i])
	                      /(2.* Delta_Y[k-1][j][i] ) ;

	Derivative_ZY[1][1] = Derivative_ZY[0][1];

	/* k+1/2*/
	Derivative_ZY[2][0] = Derivative_ZY[0][1];

	Derivative_ZY[2][1] =  0.5*(Speed_Y[k+1][j+1][i] - Speed_Y[k+1][j-1][i])
	                       /(2.* Delta_Y[k+1][j][i] ) ;
	/* i+3/2*/
	Derivative_ZY[3][0] = Derivative_ZY[0][1];

	Derivative_ZY[3][1] =  0.5*(Speed_Y[k+3][j+1][i] - Speed_Y[k+3][j-1][i])
	                      /(2.* Delta_Y[k+3][j][i] ) ;



	/* Calculation of the d/dz(du/dx) component

	  Each derivative will be interpolated 
	  in the Z-direction in order to derive the necessary quantities
	  For this reason I am using the 2D Derivative_ZX array to  compute all the interpolated
	  quantities and then combine them to get the desired result. In addition, I ll also have 
	  to interpolate in the Z-Direction, since I am using higher order schemes
	 
	  For this reason, the length Derivative_ZX array is 4 in this case
	*/

	/**** k-3/2 ****/

 	/* k-3 */
	/*\delta_1*/
	Derivative_ZX[0][0] = 2./(3.*dx)
                       	      *(Speed_X[k-3][j][i+1] - Speed_X[k-3][j][i-1] );
	/*\delta_2*/
	Derivative_ZX[0][1] = 1./(12.*dx)
	                    *(Speed_X[k-3][j][i+2] - Speed_X[k-3][j][i-2]);

	/* k */
	/*\delta_1*/
	Derivative_ZX[0][2] = 2./(3.*dx)
                       	      *(Speed_X[k][j][i+1] - Speed_X[k][j][i-1]);
	/*\delta_2*/
	Derivative_ZX[0][3] = 1./(12.*dx)
	                    *(Speed_X[k][j][i+2] - Speed_X[k-2][j][i-2]);


	/**** k-1/2 ****/	
	/*k-1*/
	/*\delta_1*/
	   Derivative_ZX[1][0] = 2./(3.*dx)           
            	      *(Speed_X[k-1][j][i+1] - Speed_X[k-1][j][i-1] );

	/*\delta_2*/
	Derivative_ZX[1][1] = 1./(12.*dx)
	                    *(Speed_X[k-1][j][i+2] - Speed_X[k-1][j][i-2]);


	/*k*/
	/*\delta_1*/	
	Derivative_ZX[1][2]=Derivative_ZX[0][2];
	/*\delta_2*/
	Derivative_ZX[1][3]=Derivative_ZX[0][3];

	/**** k+1/2 ****/
	
	/*k*/
	/*\delta_1*/
	Derivative_ZX[2][0] = Derivative_ZX[0][2];
	/*\delta_2*/
	Derivative_ZX[2][1] = Derivative_ZX[0][3];

	/*k+1*/
	/*\delta_1*/
	Derivative_ZX[2][2] = 2./(3.*dx)
                       	      *(Speed_X[k+1][j][i+1] - Speed_X[k+1][j][i-1] );
	/*\delta_2*/
	Derivative_ZX[2][3] = 1./(12.*dx)
	                    *(Speed_X[k+1][j][i+2] - Speed_X[k+1][j][i-2]);



	/**** k+3/2 ****/
	
	/*k*/
	/*\delta_1*/
	Derivative_ZX[3][0] = Derivative_ZX[0][2];
	/*\delta_2*/
	Derivative_ZX[3][1] = Derivative_ZX[0][3];

	/*k+3*/
	/*\delta_1*/
	Derivative_ZX[3][2] = 2./(3.*dx)
                       	      *(Speed_X[k+3][j][i+1] - Speed_X[k+3][j][i-1] );
	/*\delta_2*/
	Derivative_ZX[3][3] = 1./(12.*dx)
                 	    *(Speed_X[k+3][j][i+2] - Speed_X[k+3][j][i-2]);


	/* Summing the Z-component of the Viscous Term of the Z-Momentum equation*/

	/* \sigma_ij k+1/2*/
	Res_Visc_Z[2][2] = /*Viscosity(k+1/2)*/ 

	   Viscosity((Temperature[k+1][j][i]+ Temperature[k][j][i])/2. )*
	  (4./3.*(Derivative_ZZ[2][0] - Derivative_ZZ[2][1]) 
	   -2./3.*( 0.5*( (Derivative_ZY[2][1] + Derivative_ZY[2][0]) 
			+(Derivative_ZX[2][0] - Derivative_ZX[2][1])
			+(Derivative_ZX[2][2] - Derivative_ZX[2][3]) ) ));  


	/* \sigma_ij k-1/2*/
	Res_Visc_Z[2][1] = /*Viscosity(k-1/2)*/ 

	   Viscosity((Temperature[k-1][j][i]+ Temperature[k][j][i])/2 )*
	  (4./3.*(Derivative_ZZ[1][0] - Derivative_ZZ[1][1]) 
	   -2./3.*( 0.5*( (Derivative_ZY[1][1] + Derivative_ZY[1][0]) 
		       +(Derivative_ZX[1][0] - Derivative_ZX[1][1])
			+(Derivative_ZX[1][2] - Derivative_ZX[1][3]) ) ));  


	/* \sigma_ij k+3/2*/
	Res_Visc_Z[2][3]= /*Viscosity(k+3/2)*/
 
	   Viscosity((Temperature[k+1][j][i]+ Temperature[k+2][j][i])/2 )*
	  (4./3.*(Derivative_ZZ[3][0] - Derivative_ZZ[3][1]) 
	   -2./3.*( 0.5*( (Derivative_ZY[3][1] + Derivative_ZY[3][0]) 
		       +(Derivative_ZX[3][0] - Derivative_ZX[3][1])
			+(Derivative_ZX[3][2] - Derivative_ZX[3][3]) ) ));  

	/* \sigma_ij k-3/2*/
	Res_Visc_Z[2][0] =/*Viscosity(k-3/2)*/

 	   Viscosity((Temperature[k-1][j][i]+ Temperature[k-2][j][i])/2 )*
	  (4./3.*(Derivative_ZZ[0][0] - Derivative_ZZ[0][1]) 
	   -2./3.*( 0.5*( (Derivative_ZY[0][1] + Derivative_ZY[0][0]) 
		       +(Derivative_ZX[0][0] - Derivative_ZX[0][1])
		       +(Derivative_ZX[0][2] - Derivative_ZX[0][3]) ) ));  


        Res_Visc_Z_Total[2] =  9./(8.*dz)*(Res_Visc_Z[2][2] - Res_Visc_Z[2][1])
	                      -1./(24.*dz)*( Res_Visc_Z[2][3] - Res_Visc_Z[2][0]);

	    Res_Visc_Total=0;
	    for(int index=0; index<3; index++)
	      Res_Visc_Total+=Res_Visc_Z_Total[index];


	    Residual_Z[k][j][i] = -Res_Conv_Total + Res_Visc_Total+ Press_Grad;

      }
    }
  }


}
