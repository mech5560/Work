/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-05-27 14:32:34 mike_georgiou>
 *
 *
Continuity_Printer_Y.cpp -- This function will store the dv/dy
component
into an external file
*
* Written on Tuesday, 27 May 2014.
********************************************/

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <sstream>


#include"Continuity_Printer.h"

using namespace std;

void Continuity_Printer_Y(double*** velocity_y, double* dy,
                          int ldx, int ldy, int ldz,
                          int time_index, char *mike)
{

  stringstream ss;

  string format = ".dat";
  string finalName;

  ss<<mike<<time_index<<format;
  finalName = ss.str();


  ofstream outfile;
  outfile.open(finalName);


  for (int j=0; j<ldy; j++){
    for (int i=0; i<ldx; i++){

      double interpolated[2];
      for(int vj=0; vj<2; vj++)
	{
	  interpolated[vj] = 
	    Interpolation_Y(velocity_y[0][j+vj][i],dy[j+vj],
			    velocity_y[0][j+vj-1][i],dy[j+vj-1]);
	}

      double gradient_y = Derivative(interpolated[1],
				     interpolated[0],
				     2.*dy[j],1);
      
      outfile<<gradient_y<<" ";
    }

    outfile<<endl;
  }

}

