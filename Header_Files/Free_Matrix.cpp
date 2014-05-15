/*******************************************
 * Author: Michail Georgiou 
*  Last Modified: Time-stamp: <2014-04-29 18:20:37 mike_georgiou>   
*
*
Free_Matrix.cpp -- This function deallocates the memory that is being
allocated by the Matrix_Allocator function
*
* Written on Tuesday, 29 April 2014.
********************************************/

#include<cstdlib>

void Free_Matrix_New(double ***pA,int z_left, int y_left, int x_left)
{

   delete [] &pA[-z_left][-y_left][-x_left];
   delete [] &pA[-z_left][-y_left];
   delete [] &pA[-z_left];

}
