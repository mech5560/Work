/*  Last Modified Time-stamp: <2014-04-30 10:07:42 mike_georgiou> */

void Cubic_Mesh( double ***Delta_Y,
                 int ldz, int ldy, int ldx,
                 int zl, int zr,
                 int yl, int yr,
                 int xl, int xr)
{
  for (int k=0;k<ldz;k++)
    {
      for (int j=0; j<ldy;j++)
        {
          for (int i=0;i<ldx;i++)
            {
              /* Correct Version- I am just modifying it in order to evaluate the BC Function*/
              Delta_Y[k][j][i] = 0.5*(1./((double)ldy))*( (double(j+1-j)));
            }
        }
    }

  /* Extra Points necessary for the BCs */
  for (int k = 0; k < ldz; k++)
    {
      for (int j = 0; j < ldy; j++)
        {
          for (int i=1; i<=xl; i++)
            {
              /*******Left-Periodic-BC************/
                /*Delta_Y Array*/
              Delta_Y[k][j][-i] =  Delta_Y[k][j][ldx-i];
            }


          for (int i=0; i<xr; i++)
            {
              /*******Right-Periodic-BC***********/
                /*Delta_Y Array   */
              Delta_Y[k][j][ldx+i] =  Delta_Y[k][j][i];
            }

        }
    }

  /* Z-Direction BC*/
  for (int j = 0; j < ldy; j++)
    {
      for (int i = -xl; i< ldx+xr; i++)
        {
          for (int k=1; k<=zl; k++)
            {

              /******* "Back"-Periodic-BC************/
                /*Delta_Y Array*/
                Delta_Y[-k][j][i] =  Delta_Y[ldz-k][j][i];

            }

          for (int k=0; k<zr; k++)
            {

              /*******"Front"-Periodic-BC***********/
                /*Delta_Y Array   */
              Delta_Y[ldz+k][j][i] =  Delta_Y[k][j][i];
            }
        }
    }



  for (int k = -zl; k < ldz+zr; k++)
    {
      for (int i = -xl; i< ldx+xr; i++)
        {
          Delta_Y[k][-1][i] = Delta_Y[k][0][i];
          Delta_Y[k][ldy][i] = Delta_Y[k][ldy-1][i];
        }
    }

}
