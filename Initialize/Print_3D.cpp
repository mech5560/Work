/*  Last Modified Time-stamp: <2014-05-15 14:08:43 mike_georgiou> */
/* In this printing version the X-Direction is the streamwise.
   Y-Direction is the vertical
   Z-Direction is the spanwise
*/


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <sstream>


using namespace std;


void Print_3D(double ***A,
	      int ldz, int ldy, int ldx,
	      int time_index, char *mike)
{


  stringstream ss;
  string out = "";

  string format = ".dat";
  string finalName;

  ss<<out<<mike<<time_index<<format;
  finalName = ss.str();


  ofstream  outfile(finalName);



  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx;i++){

        outfile<<A[k][j][i]<<" ";
      }
      outfile<<endl;
    }
  }

}
