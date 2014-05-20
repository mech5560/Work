/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-19 16:34:52 mike_georgiou>
 *
 * BC_Velocities.cpp -- In this program I will define the
 * Boundary conditions for the three components of the velocities
 *
 *The index defines if i am having dirichlet or neuman BC
 * at the wall
 *
 *index : 1 Dirichlet
 *index : 0 Neuman
 *
 * Written on Thursday, 20 March 2014.
 **********************************************/

void BC_Velocities(double*** velocity_x,
                   double*** velocity_y,
                   double*** velocity_z,
                   int ldx, int ldy, int ldz,
                   int lx, int rx,
                   int ly, int ry,
                   int lz, int rz,
                   int index,
                   double bc_x_top, double bc_x_bottom,
                   double bc_y_top, double bc_y_bottom,
                   double bc_z_top, double bc_z_bottom,
                   double *dy)
{

  //X-Direction BC
  for (int k = 0; k < ldz; k++){
    for (int j = 0; j < ldy; j++){

      for (int i=1; i<=lx; i++)
        {
          /*******Left-Periodic-BC************/
          velocity_x[k][j][-i] =  velocity_x[k][j][ldx-i];
          velocity_y[k][j][-i] =  velocity_y[k][j][ldx-i];
          velocity_z[k][j][-i] =  velocity_z[k][j][ldx-i];
        }
      for (int i=0; i<rx; i++)
        {
          /*******Right-Periodic-BC***********/
          velocity_x[k][j][ldx+i] =  velocity_x[k][j][i];
          velocity_y[k][j][ldx+i] =  velocity_y[k][j][i];
          velocity_z[k][j][ldx+i] =  velocity_z[k][j][i];
        }

    }
  }

  // Z-Direction BC
  for (int j = 0; j < ldy; j++){
    for (int i = -lx; i< ldx+rx; i++){

      for (int k=1; k<=lz; k++)
        {
          //Back BC Periodic
          velocity_x[-k][j][i] =  velocity_x[ldz-k][j][i];
          velocity_y[-k][j][i] =  velocity_y[ldz-k][j][i];
          velocity_z[-k][j][i] =  velocity_z[ldz-k][j][i];
        }

      for (int k=0; k<rz; k++)
        {
          //Front BC Periodic
          velocity_x[ldz+k][j][i] =  velocity_x[k][j][i];
          velocity_y[ldz+k][j][i] =  velocity_y[k][j][i];
          velocity_z[ldz+k][j][i] =  velocity_z[k][j][i];
        }
    }
  }


  if (index==0) // Neuman BC = BC_Bottom
    {
      //Wall BC
      for (int k = -lz; k < ldz+rz; k++){
        for (int i = -lx; i< ldx+rx; i++){

          velocity_x[k][-ly][i] = bc_x_bottom*2.*dy[0] + velocity_x[k][0][i];
          velocity_x[k][ldy][i] = -bc_x_top*dy[ldy-1]*2. + velocity_x[k][ldy-1][i];


          velocity_y[k][-ly][i] = bc_y_bottom*2.*dy[0] + velocity_y[k][0][i];
          velocity_y[k][ldy][i] = -bc_y_top*dy[ldy-1]*2. + velocity_y[k][ldy-1][i];

          velocity_z[k][-ly][i] = bc_z_bottom*2.*dy[0] + velocity_z[k][0][i];
          velocity_z[k][ldy][i] = -bc_z_top*dy[ldy-1]*2. + velocity_z[k][ldy-1][i];


        }
      }
    }

  if (index==1) // Dirichlet BC = BC_Bottom
    {

      for (int k = -lz; k < ldz+rz; k++){
        for (int i = -lx; i< ldx+rx; i++)
          {

            velocity_x[k][-ly][i] = 2.*bc_x_bottom - velocity_x[k][0][i];
            velocity_x[k][ldy][i] = 2.*bc_x_top - velocity_x[k][ldy-1][i];

            velocity_y[k][-ly][i] = 2.*bc_y_bottom - velocity_y[k][0][i];
            velocity_y[k][ldy][i] = 2.*bc_y_top - velocity_y[k][ldy-1][i];

            velocity_z[k][-ly][i] = 2.*bc_z_bottom - velocity_z[k][0][i];
            velocity_z[k][ldy][i] = 2.*bc_z_top - velocity_z[k][ldy-1][i];

          }
      }
    }
}
