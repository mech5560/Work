/*  Last Modified Time-stamp: <2014-04-30 09:14:48 mike_georgiou> */
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;
void Print_3D_Single(double *A,
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


  // for (int k=0; k<ldz; k++)
  //   {
  for (int j=0; j<ldy; j++)
    {
      for (int i=0; i<ldx; i++)
        {
          outfile<<A[(0)*ldy*ldx + ldy*j +i+1]<<" ";
        }
      outfile<<endl;
    }
  outfile<<endl;

  //    }
}
