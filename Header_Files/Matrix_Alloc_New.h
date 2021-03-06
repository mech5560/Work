/*  Last Modified Time-stamp: <2014-03-19 17:22:24 mike_georgiou> */ 
#ifndef MATRIX_ALLOC_NEW_H
#define MATRIX_ALLOC_NEW_H

inline double *** Matrix_Alloc_New(int ldz, int ldy, int ldx, 
			     int z_left,int z_right,
			     int y_left,int y_right,
			     int x_left,int x_right)
{

  int z_tot=ldz+z_left+z_right;
  int y_tot=ldy+y_left+y_right;
  int x_tot=ldx+x_left+x_right;

  double ***pA;
  
  pA= new double**[z_tot]; pA+=z_left;
  pA[-z_left]= new double* [z_tot*y_tot] ; pA[-z_left]+=y_left;
  pA[-z_left][-y_left]= new double [z_tot*y_tot*x_tot] ; pA[-z_left][-y_left] +=x_left;


  
  // Pointer Assignment
  for (int j = -z_left; j < ldz+z_right; j++) 
    {
      if (j !=-z_left)
  	{
  	  pA[j] = pA[j-1] + y_tot;
  	  pA[j][-y_left] = pA[j-1][-y_left] + y_tot*x_tot;
  	}
      for (int i=-y_left+1; i<ldy+y_right; i++)
  	{
  	  pA[j][i] = pA[j][i-1] + x_tot;
  	}
      
    }

  return pA;


}
#endif

#ifndef MATRIX_FREE_NEW
#define MATRIX_FREE_NEW

inline void Free_Matrix_New(double ***pA,int z_left, int y_left, int x_left)
{

   delete [] &pA[-z_left][-y_left][-x_left];
   delete [] &pA[-z_left][-y_left];
   delete [] &pA[-z_left];

}

#endif
