/*******************************************
* Author: Michail Georgiou
*  Last Modified: Time-stamp: <2014-03-25 10:01:48 mike_georgiou>
*
* Boundary_Conditions.cpp -- In this program I will define  the boundary conditions for my problem.
* Firstly, I ll start with the easy ones (Periodic) and then I ll define the BC at the wall
   *
   * The velocity BC at the wall is defined by the non-slip condition.
   * For the density I can define either Dirichlet=0 or obtain the Rho_{wall-1} by the constitutive equatiion.
   * As for the intermediate velocities i ll follow christos suggestions
   *
   * Written on Thursday, 20 March 2014.
   ********************************************/


#include "Struct.h"
#include "./../../Header_Files/Data.h"
 
  void BC(Ar *Arr, int ldz, int ldy, int ldx,
          int lz, int rz, int ly, int ry, int lx, int rx, int index)
  {


    /* X- Direction BC */
    for (int k = 0; k < ldz; k++)
      {
        for (int j = 0; j < ldy; j++)
          {
            for (int i=1; i<=lx; i++)
              {
                /*******Left-Periodic-BC************/

		/*Delta_Y Array*/
		Arr->Delta_Y[k][j][-i] =  Arr->Delta_Y[k][j][ldx-i];

                /* Velocities Array */
                Arr->Speed_X[k][j][-i] =  Arr->Speed_X[k][j][ldx-i];
                Arr->Speed_Y[k][j][-i] =  Arr->Speed_Y[k][j][ldx-i];
                Arr->Speed_Z[k][j][-i] =  Arr->Speed_Z[k][j][ldx-i];

                Arr->Speed_New_X[k][j][-i] =  Arr->Speed_New_X[k][j][ldx-i];
                Arr->Speed_New_Y[k][j][-i] =  Arr->Speed_New_Y[k][j][ldx-i];
                Arr->Speed_New_Z[k][j][-i] =  Arr->Speed_New_Z[k][j][ldx-i];

                /*Densities*/
                Arr->Rho[k][j][-i] = Arr->Rho[k][j][ldx-i];
                Arr->Rho_Nm1[k][j][-i] = Arr->Rho_Nm1[k][j][ldx-i];
                Arr->Rho_New[k][j][-i] = Arr->Rho_New[k][j][ldx-i];

                /* Pressure*/
                Arr->Pressure[k][j][-i] = Arr->Pressure[k][j][ldx-i];
                Arr->Pressure_New[k][j][-i] = Arr->Pressure_New[k][j][ldx-i];

                /*Temperature*/
                Arr->Temperature[k][j][-i] = Arr->Temperature[k][j][ldx-i];
		//Arr->Pressure_New[k][j][-i] = Arr->Pressure_New[k][j][ldx-i];






              }

            for (int i=0; i<rx; i++)
              {

                /*******Right-Periodic-BC***********/

                /*Delta_Y Array   */
                Arr->Delta_Y[k][j][ldx+i] =  Arr->Delta_Y[k][j][i];

                /* Velocities Array */
                Arr->Speed_X[k][j][ldx+i] =  Arr->Speed_X[k][j][i];
                Arr->Speed_Y[k][j][ldx+i] =  Arr->Speed_Y[k][j][i];
                Arr->Speed_Z[k][j][ldx+i] =  Arr->Speed_Z[k][j][i];

                Arr->Speed_New_X[k][j][ldx+i] =  Arr->Speed_New_X[k][j][i];
                Arr->Speed_New_Y[k][j][ldx+i] =  Arr->Speed_New_Y[k][j][i];
                Arr->Speed_New_Z[k][j][ldx+i] =  Arr->Speed_New_Z[k][j][i];

                /*Densities*/
                Arr->Rho[k][j][ldx+i] = Arr->Rho[k][j][i];
                Arr->Rho_Nm1[k][j][ldx+i] = Arr->Rho_Nm1[k][j][i];
                Arr->Rho_New[k][j][ldx+i] = Arr->Rho_New[k][j][i];

                /* Pressure*/
                Arr->Pressure[k][j][ldx+i] = Arr->Pressure[k][j][i];
                Arr->Pressure_New[k][j][ldx+i] = Arr->Pressure_New[k][j][i];
		
		/*Temperature */
		Arr->Temperature[k][j][ldx+i] = Arr->Temperature[k][j][i];

              }

          }
      }


    /* Z-Direction BC*/
    for (int j = 0; j < ldy; j++)
      {
        for (int i = -lx; i< ldx+rx; i++)
          {
            for (int k=1; k<=lz; k++)
              {
    
		/******* "Back"-Periodic-BC************/
		  
    		/*Delta_Y Array*/
    		Arr->Delta_Y[-k][j][i] =  Arr->Delta_Y[ldz-k][j][i];

                /* Velocities Array */
                Arr->Speed_X[-k][j][i] =  Arr->Speed_X[ldz-k][j][i];
                Arr->Speed_Y[-k][j][i] =  Arr->Speed_Y[ldz-k][j][i];
                Arr->Speed_Z[-k][j][i] =  Arr->Speed_Z[ldz-k][j][i];

                Arr->Speed_New_X[-k][j][i] =  Arr->Speed_New_X[ldz-k][j][i];
                Arr->Speed_New_Y[-k][j][i] =  Arr->Speed_New_Y[ldz-k][j][i];
                Arr->Speed_New_Z[-k][j][i] =  Arr->Speed_New_Z[ldz-k][j][i];

                /*Densities*/
                Arr->Rho[-k][j][i] = Arr->Rho[ldz-k][j][i];
                Arr->Rho_Nm1[-k][j][i] = Arr->Rho_Nm1[ldz-k][j][i];
                Arr->Rho_New[-k][j][i] = Arr->Rho_New[ldz-k][j][i];

                /* Pressure*/
                Arr->Pressure[-k][j][i] = Arr->Pressure[ldz-k][j][i];
                Arr->Pressure_New[-k][j][i] = Arr->Pressure_New[ldz-k][j][i];

		/*Temperature*/
		Arr->Temperature[-k][j][i] = Arr-> Temperature[ldz-k][j][i];

              }

            for (int k=0; k<rz; k++)
              {

                /*******"Front"-Periodic-BC***********/

                /*Delta_Y Array   */
                Arr->Delta_Y[ldz+k][j][i] =  Arr->Delta_Y[k][j][i];

                /* Velocities Array */
                Arr->Speed_X[ldz+k][j][i] =  Arr->Speed_X[k][j][i];
                Arr->Speed_Y[ldz+k][j][i] =  Arr->Speed_Y[k][j][i];
                Arr->Speed_Z[ldz+k][j][i] =  Arr->Speed_Z[k][j][i];

                Arr->Speed_New_X[ldz+k][j][i] =  Arr->Speed_New_X[k][j][i];
                Arr->Speed_New_Y[ldz+k][j][i] =  Arr->Speed_New_Y[k][j][i];
                Arr->Speed_New_Z[ldz+k][j][i] =  Arr->Speed_New_Z[k][j][i];

                /*Densities*/
                Arr->Rho[ldz+k][j][i] = Arr->Rho[k][j][i];
                Arr->Rho_Nm1[ldz+k][j][i] = Arr->Rho_Nm1[k][j][i];
                Arr->Rho_New[ldz+k][j][i] = Arr->Rho_New[k][j][i];

                /* Pressure*/
                Arr->Pressure[ldz+k][j][i] = Arr->Pressure[k][j][i];
                Arr->Pressure_New[ldz+k][j][i] = Arr->Pressure_New[k][j][i];
		
		/*Temperature*/
		Arr->Temperature[ldz+k][j][i] = Arr-> Temperature[k][j][i];


              }
          }
      }

    /* Wall BC for all the quantities */
			 

    for (int k = -lz; k < ldz+rz; k++)
      {
        for (int i = -lx; i< ldx+rx; i++)
          {
	  
	    /* Velocities - No slip Condition */
	 
	    /*Bottom Wall*/
	    Arr->Speed_X[k][-ly][i] = -Arr->Speed_X[k][0][i];
	    Arr->Speed_Y[k][-ly][i] = -Arr->Speed_Y[k][0][i];
	    Arr->Speed_Z[k][-ly][i] = -Arr->Speed_Z[k][0][i];
	
	    /*Top Wall*/
	    Arr->Speed_X[k][ldy][i] = -Arr->Speed_X[k][ldy-1][i];
	    Arr->Speed_Y[k][ldy][i] = -Arr->Speed_Y[k][ldy-1][i];
	    Arr->Speed_Z[k][ldy][i] = -Arr->Speed_Z[k][ldy-1][i];



	    /* Temperatures - for now i ll put a dirichlet BC 
	    and i ll have uniform temperature everywhere*/ 

	    /*Bottom Wall*/
	    Arr->Temperature[k][-ly][i] = 2.*Temp_Bottom - Arr->Temperature[k][0][i];

	    /*Top Wall*/
	    Arr->Temperature[k][ldy][i] = 2.*Temp_Top - Arr->Temperature[k][ldy-1][i];


	    /*Densities - For now I ll set Neuman=0 for boundary condition*/

	    /* Bottom Wall*/
	    Arr->Rho_Nm1[k][-ly][i] = 2.*Arr->Delta_Y[k][0][i] * Rho_Grad_Bottom + Arr->Rho_Nm1[k][0][i];
	    Arr->Rho[k][-ly][i] = 2.*Arr->Delta_Y[k][0][i] * Rho_Grad_Bottom + Arr->Rho[k][0][i];

	    /* Top  Wall*/
	    Arr->Rho_Nm1[k][ldy][i] = 2.*Arr->Delta_Y[k][ldy-1][i] * Rho_Grad_Top + Arr->Rho_Nm1[k][ldy-1][i];
	    Arr->Rho[k][ldy][i] = 2.*Arr->Delta_Y[k][ldy-1][i] * Rho_Grad_Top + Arr->Rho[k][ldy-1][i];

	  }
      }





    
}