/*  Last Modified Time-stamp: <2014-05-16 09:44:14 mike_georgiou> */

void Cubic_Mesh( double *dy, double *y,
		 double length_y,
                 int ldy,
                 int yl, int yr)
{
  double dy_local= length_y/ldy;

  for (int j=0; j<ldy;j++)
    {
      // Correct Version- I am just modifying it in order
      //         to evaluate the BC Function
      dy[j] = dy_local/2.;
    }

  double y_local=-1.;
  for (int j=0; j<ldy;j++)
    {
      y_local += dy[j];
      y[j] = y_local;
      y_local += dy[j];
    }



  //Boundary Conditions
  dy[-1] = dy[0];
  dy[ldy] = dy[ldy-1];

  y[-1] = y[0];
  y[ldy] = y[ldy-1];


}
