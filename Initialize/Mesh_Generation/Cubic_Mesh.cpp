/*  Last Modified Time-stamp: <2014-05-01 16:05:56 mike_georgiou> */

void Cubic_Mesh( double *dy,
                 int ldy,
                 int yl, int yr)
{
      for (int j=0; j<ldy;j++)
        {
					// Correct Version- I am just modifying it in order 
					// 	 to evaluate the BC Function
					dy[j] = 0.5*(1./((double)ldy))*( (double(j+1-j)));
				}

			//Boundary Conditions			
          dy[-1] = dy[0];
          dy[ldy] = dy[ldy-1];

}
