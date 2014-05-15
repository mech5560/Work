/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-01 17:42:20 mike_georgiou>
 *
 * Boundary_Conditions.cpp -- In this program I will define  the boundary
 * conditions  for my problem.
 * Firstly, I ll start with the easy ones (Periodic) and then I ll define the
 * BC at the wall
 *
 * The velocity BC at the wall is defined by the non-slip condition.
 * For the density I can define either Dirichlet=0 or obtain the Rho_{wall-1} by 
the constitutive equatiion.
 * As for the intermediate velocities i ll follow christos suggestions
 *
 * Written on Thursday, 20 March 2014.
 ********************************************/


#include "Struct.h"
#include "./../../Header_Files/Data.h"

void BC_Predictor(Ar *Arr,
                  int ldz, int ldy, int ldx,
                  int lz, int rz,
                  int ly, int ry,
                  int lx, int rx,
                  int index)
{


  /* X- Direction BC */
  for (int k = 0; k < ldz; k++){
    for (int j = 0; j < ldy; j++){
      for (int i=1; i<=lx; i++){

        /*******Left-Periodic-BC************/


        /* Velocities Array */
        Arr->velocity_x[k][j][-i] =  Arr->velocity_x[k][j][ldx-i];
        Arr->velocity_y[k][j][-i] =  Arr->velocity_y[k][j][ldx-i];
        Arr->velocity_z[k][j][-i] =  Arr->velocity_z[k][j][ldx-i];

      }

      for (int i=0; i<rx; i++)
        {

          /*******Right-Periodic-BC***********/

          /* Velocities Array */
          Arr->velocity_x[k][j][ldx+i] =  Arr->velocity_x[k][j][i];
          Arr->velocity_y[k][j][ldx+i] =  Arr->velocity_y[k][j][i];
          Arr->velocity_z[k][j][ldx+i] =  Arr->velocity_z[k][j][i];

        }

    }
  }


  /* z-Direction BC*/
  for (int j = 0; j < ldy; j++){
    for (int i = -lx; i< ldx+rx; i++){
      for (int k=1; k<=lz; k++){

        /******* "Back"-Periodic-BC************/

        /* Velocities Array */
        Arr->velocity_x[-k][j][i] =  Arr->velocity_x[ldz-k][j][i];
        Arr->velocity_y[-k][j][i] =  Arr->velocity_y[ldz-k][j][i];
        Arr->velocity_z[-k][j][i] =  Arr->velocity_z[ldz-k][j][i];

      }

      for (int k=0; k<rz; k++)
        {


          /* Velocities Array */
          Arr->velocity_x[ldz+k][j][i] =  Arr->velocity_x[k][j][i];
          Arr->velocity_y[ldz+k][j][i] =  Arr->velocity_y[k][j][i];
          Arr->velocity_z[ldz+k][j][i] =  Arr->velocity_z[k][j][i];


        }
    }
  }

  /* Wall BC for all the quantities */


  for (int k = -lz; k < ldz+rz; k++){
    for (int i = -lx; i< ldx+rx; i++){

      /* Velocities - No slip Condition */

      /*Bottom Wall*/
      Arr->velocity_x[k][-ly][i] = -Arr->velocity_x[k][0][i];
      Arr->velocity_y[k][-ly][i] = -Arr->velocity_y[k][0][i];
      Arr->velocity_z[k][-ly][i] = -Arr->velocity_z[k][0][i];

      /*Top Wall*/
      Arr->velocity_x[k][ldy][i] = -Arr->velocity_x[k][ldy-1][i];
      Arr->velocity_y[k][ldy][i] = -Arr->velocity_y[k][ldy-1][i];
      Arr->velocity_z[k][ldy][i] = -Arr->velocity_z[k][ldy-1][i];

    }
  }






}
