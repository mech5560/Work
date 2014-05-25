/*******************************************
 * Author: Michail Georgiou
 * Last Modified: Time-stamp: <2014-05-22 15:47:10 mike_georgiou>
 *
 * BC_Single.cpp -- In this program I will define the
 * Boundary conditions for the Density and the Pressure.
 *
 *The index defines if i am having dirichlet or neuman BC
 * at the wall
 *
 *index : 1 Dirichlet
 *index : 2 Neuman
 *
 * Written on Thursday, 20 March 2014.
 **********************************************/


//#include "./../../Header_Files/Data.h"

void BC_Single(double ***data,
               int ldx, int ldy, int ldz,
               int lx, int rx,
               int ly, int ry,
               int lz, int rz,
               int index,
               double bc_top, double bc_bottom,
               double *dy)
{

  //X-Direction BC
  for (int k = 0; k < ldz; k++){
    for (int j = 0; j < ldy; j++){

      for (int i=1; i<=lx; i++)
        {
          /*******Left-Periodic-BC************/
          data[k][j][-i] = data[k][j][ldx-i];
        }
      for (int i=0; i<rx; i++)
        {
          /*******Right-Periodic-BC***********/
          data[k][j][ldx+i] = data[k][j][i];
        }

    }
  }

  // Z-Direction BC
  for (int j = 0; j < ldy; j++){
    for (int i = -lx; i< ldx+rx; i++){

      for (int k=1; k<=lz; k++)
        {
          //Back BC Periodic
          data[-k][j][i] = data[ldz-k][j][i];
        }

      for (int k=0; k<rz; k++)
        {
          //Front BC Periodic
          data[ldz+k][j][i] = data[k][j][i];
        }
    }
  }


  if (index==2) // Neuman BC = BC_Bottom
    {
      //Wall BC
      for (int k = -lz; k < ldz+rz; k++){
        for (int i = -lx; i< ldx+rx; i++){
            ///NEED TO FIX THAT ERROR WITH DY
          data[k][-ly][i] = bc_bottom*2.*dy[0] + data[k][0][i];
          data[k][ldy][i] = -bc_top*dy[ldy-1]*2. + data[k][ldy-1][i];
        }
      }
    }

  if (index==1) // Dirichlet BC = BC_Bottom
    {

      for (int k = -lz; k < ldz+rz; k++){
        for (int i = -lx; i< ldx+rx; i++)
          {

            data[k][-ly][i] = 2.*bc_bottom - data[k][0][i];
            data[k][ldy][i] = 2.*bc_top - data[k][ldy-1][i];
          }
      }
    }
}
