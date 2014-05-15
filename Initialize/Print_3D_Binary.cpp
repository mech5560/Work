/*  Last Modified Time-stamp: <2014-04-30 11:21:12 mike_georgiou> */
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;

void Print_3D_Binary(double ***A,
                     int ldz, int ldy, int ldx,
                     int time_index, char *mike)
{


  stringstream ss;
  string out = "";

  string format = ".bin";
  string finalName;

  ss<<out<<mike<<time_index<<format;
  finalName = ss.str();

  ofstream  outfile;
  outfile.open(finalName,ios::out | ios::binary);


  for (int j=0; j<ldy; j++)
    {
      for (int i=0; i<ldx; i++)
        {
          outfile.write( (char*)&A[4][j][i],sizeof(double));
        }
    }

}
