/*  Last Modified Time-stamp: <2014-03-26 17:29:46 mike_georgiou> */ 
/* In this printing version the X-Direction is the streamwise. 
                                Y-Direction is the vertical
                                Z-Direction is the spanwise
*/

#include <iostream>
#include <stdio.h>
using namespace std;

void Print_3D(double ***A, int ldz, int ldy, int ldx,
	      int zl,int  zr,
	      int yl, int yr, 
	      int xl, int xr)
{

  int input;
  cout<< "Pause? (1: for yes, 0: for no)\n";
    cin>>input;

    for (int k=-zl; k<ldz+zr; k++)
      {
	for (int j=-yl; j<ldy+yr; j++)
	  {
	    for (int i=-xl; i<ldx+xr; i++)
	      {
		cout<<A[k][j][i]<<" ";
	      }
	    cout<<endl;
	  }
	cout<<endl;
      if (input==1)
	getchar();
    }
}
