/*******************************************
 * Author: Michail Georgiou 
*  Last Modified: Time-stamp: <2014-04-30 17:35:00 mike_georgiou>   
*
*
Free_Matrix.cpp -- This function deallocates the memory that is being
allocated by the Matrix_Allocator function
*
* Written on Tuesday, 29 April 2014.
********************************************/

 void Free_Matrix(double ***pA,
									int x_left, int y_left, int z_left)
{

   delete [] &pA[-z_left][-y_left][-x_left];
   delete [] &pA[-z_left][-y_left];
   delete [] &pA[-z_left];

}
