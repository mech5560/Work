/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-15 13:05:51 mike_georgiou>
 *
 *
Print_2D_Data.cpp -- This program prints writes two columns of data that are
necessary for gnuplot to visualize the results of a non-uniform grid.
*
* Written on Friday,  2 May 2014.
********************************************/
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;

void Print_2D_Data(double*** data, double* dy,
                   int ldx, int ldy, int ldz,
                   int time_index, char* name)
{

  stringstream ss;
  string out = "";

  string format = ".dat";
  string finalName;

  ss<<out<<name<<time_index<<format;
  finalName = ss.str();


  ofstream  outfile(finalName);

  outfile<<"#  X     Y"<<endl;
  // for (int k=0; k<ldz; k++){
  double y_coord=-1.;
  for (int j=0; j<ldy; j++)
    {
    y_coord+=dy[j];
    outfile<<y_coord<<" "<<data[0][j][0]<<endl;
    y_coord+=dy[j];
    }

  outfile<<endl;
}
