/*  Last Modified Time-stamp: <2014-04-29 18:20:42 mike_georgiou> */ 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>


#ifndef Allocate_Memory
#define Allocate_Memory

double *** Matrix_Alloc(int ldz, int ldy, int ldx)
{
  double*** Speed_X;
  
  Speed_X = new double**[ldz];
  Speed_X[0] = new double*[ldz * ldy];
  Speed_X[0][0] = new double[ldz * ldy * ldx ];

  for (int i = 0; i < ldz; i++)
    {
      if (i > 0)
	{
	  Speed_X[i] = Speed_X[i-1] + ldy;
	  Speed_X[i][0] = Speed_X[i-1][0] + (ldy*ldx);
	}
      for (int j = 1; j < ldy; j++)
	{
	  Speed_X[i][j] = Speed_X[i][j-1] + ldx;
	}
    }

  return Speed_X;
}
#endif

#ifndef Free_Memory
#define Free_Memory
 void Free_Matrix(double ***Speed_X, int ldz, int ldy, int ldx)
{
  delete [] &Speed_X[0][0][0];
  delete [] &Speed_X[0][0];
  delete [] &Speed_X[0];
 
}
#endif




