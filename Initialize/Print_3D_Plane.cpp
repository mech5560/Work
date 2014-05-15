/*  Last Modified Time-stamp: <2014-03-27 12:45:13 mike_georgiou> */ 
/* In this printing version the X-Direction is the streamwise. 
                                Y-Direction is the vertical
                                Z-Direction is the spanwise
*/

#include <iostream>
#include <stdio.h>
using namespace std;

void Print_3D_Plane(double ***A, int ldz, int ldy, int ldx,
	      int zl,int  zr,
	      int yl, int yr, 
	      int xl, int xr)
{



	for (int j=-yl; j<ldy+yr; j++)
	  {
	    for (int i=-xl; i<ldx+xr; i++)
	      {
		cout<<A[ldz-1][j][i]<<" ";
	      }
	    cout<<endl;
	  }

	getchar();
}
