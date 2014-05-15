/*  Last Modified Time-stamp: <2014-04-24 16:15:30 mike_georgiou> */ 
/* Definition of macros that will be used on the code */

//#include "Data.h"
#include <iostream>
#include <cstdlib>
#include <math.h>
using namespace std;

//1 Dynamic Memory Allocation Controller
#ifndef Alloc_Check 
#define  Alloc_Check(host_A_vec)					\
  if ((host_A_vec) == 0){						\
    cout << "Error: memory could not be allocated in this Array"<<endl;	\
    exit(1); }
#endif


#ifndef A
// Macros that "Converts" 1D array into 3D
#define A(z,y,x) ((z)*ldx*ldy + (y)*ldx + (x))	
#endif


#ifndef A_sp
// Macros that "Converts" 1D array into 3D
#define A_sp(z,y,x) ((z)*ldx*(ldy+2) + (y)*ldx + (x))	
#endif



/* I ll use these function to avoid multiple if statements close to the 
   boundaries in the X and Z Directions
*/

#ifndef A_Bound
#define A_Bound 

inline  int A_Boundary_ldz(int k,int ldz)
{
  if(k>=ldz){k=k-(ldz);}
  return k;
}

inline  int A_Boundary_0(int k,int ldz)
{
  if(k<0){ k= ldz+k;}
  return k;
}

#endif


#ifndef Periodic_Check
#define Periodic_Check
inline void Periodic_Bound_Check(int i, int *i_ref, int ldx)
{

	
	if (i<3)
	  {
	    i_ref[0]=A_Boundary_0(i_ref[0],ldx);
	    i_ref[1]=A_Boundary_0(i_ref[1],ldx);
	    i_ref[2]=A_Boundary_0(i_ref[2],ldx);	 
	  }

	if (i>ldx-4)
	  {
	    i_ref[3]=A_Boundary_ldz(i_ref[3],ldx);
	    i_ref[4]=A_Boundary_ldz(i_ref[4],ldx);
	    i_ref[5]=A_Boundary_ldz(i_ref[5],ldx);	 
	  }


}
#endif



/**********************************************************
 The next two formulas must be non-dimensionalized with the 
 reference values

 Relations found by "The Shock Absorber Handbook (Properties of Water)"

**********************************************************/

/* //4 Calculation of Thermal Conductivity */
/* #ifndef Conductivity */
/* #define Conductivity(T)				\ */
/*   0.886025 + 0.001756/0.644*(T) - 0.00000646/0.644*(T)*(T) */
/* #endif */


/* //4 Calculation of Heat Capacity */
/* #ifndef Heat_Capacity */
/* #define Heat_Capacity(T)			\ */
/*   1.005495 -1.31/4186.0*(T) +0.014/4186.*(T)*(T) */
/* #endif */


// 5 Calculation of Density
// Alternatively, I can use a relationship with constant volumetric expansion
#ifndef Density
#define Density(T) \
1.034615 -0.155/988.*(T) - 0.0026580/988. * (T)*(T)
#endif

//6 Calculation of Viscosity
#ifndef Viscosity
#define Viscosity(T)\
1. 
 //  pow(10, (-2.75 - 0.0141*(T) +91.9*1e-6*(T)*(T) -311*1e-9 *(T)*(T)*(T))/(0.547e-3) )
#endif
