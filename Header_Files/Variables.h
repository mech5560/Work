/*  Last Modified Time-stamp: <2014-04-30 16:59:29 mike_georgiou> */ 
/* In this header file I will define the variables that
   I will use in this project
*/
//#include"Struct.h"

#ifndef Variables
#define Variables


int 
  ldx,ldy,ldz,
/*Nze is defined in order to construct the vectors for the sparse solver*/
  nze, dim_a, flag;

int left_z,right_z,
	left_y,right_y,
	left_x,right_x;



/* Structure defined in the Struct.h file*/
Ar Arr;


double *s_a, *precond_a, *result, *rhs,Coefficients_X[4],Coefficients_Z[4], ***Temp;
int *ij_a;
double dx, dz;
#endif
