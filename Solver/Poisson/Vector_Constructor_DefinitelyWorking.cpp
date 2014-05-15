
struct Constructor
{

  int Position[14][2];
  double Value [14];
  double Diagonal;
  double dy[3];
};

struct Constructor_Y
{

  int Position_Y[13][2];
  double Value_Y[13];
  double Diagonal;
  double dy_special[2];
};

#include "./../../Header_Files/Macros.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

void bubbleSort( int size,int array[][2]);

/*
  When I ll use this function in my code I ll have to modify it in order to accept a triple pointer
  instead of a single one for the Delta_Y array
*/

void Vector_Constructor(double *s_A, double *Precond_A, int *ij_A,
                        double Coefficients_Z[4], double Coefficients_X[4], double ***Delta_Y,
                        int ldz, int ldy, int ldx, int Nze)
{


  struct Constructor General_Case;
  struct Constructor_Y Special_Case;

  int k_ref[6],i_ref[6];
  int Element_Counter;


  /* There is a chance that I have an error here  DOUBLECHECK LATER*/
  Element_Counter = ldz*ldy*ldx +2;

  ij_A[Element_Counter-1] = Nze+2;
  s_A[Element_Counter-1] = 0.;


  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){



        /* passing the stencil indices to the k_ref and i_ref arrays
           in order to treat the periodic B.C in the X and Z direction periodically*/


        for (int index=0; index<3; index++)
          {
            i_ref[index] = i-(3-index);
            i_ref[index+3] = i+(1+index);

            k_ref[index] = k-(3-index);
            k_ref[index+3] = k+(1+index);

          }

        /*Checking if I am close to the Boundaries in the X and Z directions*/
        Periodic_Bound_Check(k, k_ref, ldz);
        Periodic_Bound_Check(i, i_ref, ldx);


        if (j>0 && j<ldy-1)
          {


            /* Calculating the second derivatives coefficients in the non-uniform direction*/

	      /*j-1 Coefficient*/
	    General_Case.dy[0]= 1./(2.*Delta_Y[k][j][i]
	    *(Delta_Y[k][j+1][i]+Delta_Y[k][j][i]));
				     
	    General_Case.dy[1] =
	      -(Delta_Y[k][j+1][i]+2.*Delta_Y[k][j][i]+Delta_Y[k][j-1][i])/(
	      (2.*Delta_Y[k][j][i])*(Delta_Y[k][j+1][i]+Delta_Y[k][j][i])
	      *(Delta_Y[k][j-1][i]+Delta_Y[k][j][i]));
				    
            /*j+1 Coefficient*/ 
            General_Case.dy[2] = 1./(2.*Delta_Y[k][j][i]
	    *(Delta_Y[k][j-1][i]+Delta_Y[k][j][i]));

	    // General_Case.dy[0]=1./(dy*dy);
	    // General_Case.dy[1]=-2./(dy*dy);
	    // General_Case.dy[2]=1./(dy*dy);

            /*Defining the nze position in the A (square) matrix*/
            for (int index=0; index<3; index++)
              {

                /*k-3:k-1*/
                General_Case.Position[index][0] = A(k_ref[index],j,i) +1;
                General_Case.Position[index][1] = index;
                General_Case.Value[index] = Coefficients_Z[index];

                /*j-1 outside loop*/

                /*i-3:i-1*/
                General_Case.Position[index+4][0] = A(k,j,i_ref[index]) +1;
                General_Case.Position[index+4][1] = index+4;
                General_Case.Value[index+4] = Coefficients_X[index];

                /* Diagonal Element*/

                /*i+1:i+3*/
                General_Case.Position[index+7][0] = A(k,j,i_ref[index+3]) +1;
                General_Case.Position[index+7][1] = index+7;
                General_Case.Value[index+7] = Coefficients_X[2-index];

                /*j+1 outside loop*/

                /*k+1:k+3*/
                General_Case.Position[index+11][0] = A(k_ref[index+3],j,i) +1;
                General_Case.Position[index+11][1] = index+11;
                General_Case.Value[index+11] = Coefficients_Z[2-index];

              }


            /* j-1*/
            General_Case.Position[3][0] = A(k,j-1,i) +1;
            General_Case.Position[3][1] = 3;
            General_Case.Value[3] = General_Case.dy[0];


            /*j+1*/
            General_Case.Position[10][0] = A(k,j+1,i) +1;
            General_Case.Position[10][1] = 10;
            General_Case.Value[10] = General_Case.dy[2];

            General_Case.Diagonal = Coefficients_X[3] + Coefficients_Z[3]
              +General_Case.dy[1];
	    



            /* Calling the sorting function */
            bubbleSort(14,General_Case.Position);


            /* Filling the first 1:Dim_A elements of the s_A vector with the
               diagonal elements */
            s_A[A(k,j,i)+1] = General_Case.Diagonal;

            /* Filling the Preconditioner vector with the diagonal Values of A*/
            Precond_A[A(k,j,i) + 1] =  General_Case.Diagonal;

            /* Filling the first 1:Dim_A elements of the ij_A vector with the
               position of the first off diagonal elements of each row*/
            ij_A [A(k,j,i)+1] = Element_Counter;


            /*
              Filling the Dim_A+2 : NZE  of the s_A and ij_A matrices with
              the off diagonal elements and their corresponding column indices
            */
            for (int index=0; index<14; index++)
              {
                /* Value of the Rest off-Diagonal Elements */
                s_A[Element_Counter +index] =
                  General_Case.Value[General_Case.Position[index][1]];

                /* Column index of the corresponding elements*/
                ij_A [Element_Counter + index] = General_Case.Position[index][0];
              }

            Element_Counter += 14;

          }

         if(j==0)
          {
            Special_Case.dy_special[0] = -1./
               1./(2.*Delta_Y[k][j][i]
		   *(Delta_Y[k][j+1][i]+Delta_Y[k][j][i]));
           
            Special_Case.dy_special[1] = 1./
               1./(2.*Delta_Y[k][j][i]
		   *(Delta_Y[k][j+1][i]+Delta_Y[k][j][i]));
	    


            /*testing Purpuses */

	    // Special_Case.dy_special[0]=-1./(dy*dy);
            // Special_Case.dy_special[1]=1./(dy*dy);

	       
	   

            for (int index=0; index<3; index++)
              {

                /*k-3:k-1*/
                Special_Case.Position_Y[index][0] = A(k_ref[index],j,i) +1;
                Special_Case.Position_Y[index][1] = index;
                Special_Case.Value_Y[index] = Coefficients_Z[index];

                /*i-3:i-1*/
                Special_Case.Position_Y[index+3][0] = A(k,j,i_ref[index]) +1;
                Special_Case.Position_Y[index+3][1] = index+3;
                Special_Case.Value_Y[index+3] = Coefficients_X[index];

                /* Diagonal Element*/

                /*i+1:i+3*/
                Special_Case.Position_Y[index+6][0] = A(k,j,i_ref[index+3]) +1;
                Special_Case.Position_Y[index+6][1] = index+6;
                Special_Case.Value_Y[index+6] = Coefficients_X[2-index];

                /*j+1 outside loop*/

                /*k+1:k+3*/
                Special_Case.Position_Y[index+10][0] = A(k_ref[index+3],j,i) +1;
                Special_Case.Position_Y[index+10][1] = index+10;
                Special_Case.Value_Y[index+10] = Coefficients_Z[2-index];

              }

            /*j+1*/
            Special_Case.Position_Y[9][0] = A(k,j+1,i) +1;
            Special_Case.Position_Y[9][1] = 9;
            Special_Case.Value_Y[9] = Special_Case.dy_special[1];


            /* Diagonal Element */

            Special_Case.Diagonal = Coefficients_X[3] + Coefficients_Z[3]
              +Special_Case.dy_special[0];



            /* Calling the Sorting Function*/
            bubbleSort(13,Special_Case.Position_Y);



            /* Filling the first 1:Element_Counter elements of the s_A vector with the
               diagonal elements */
            s_A[A(k,j,i)+1] = Special_Case.Diagonal;

            /* Filling the Preconditioner vector with the diagonal Values of A*/
            Precond_A[A(k,j,i) + 1] =  Special_Case.Diagonal;


            /* Filling the first 1:Element_Counter elements of the ij_A vector with the
               position of the first off diagonal elements of each row*/
            ij_A [A(k,j,i)+1] = Element_Counter;


            /*
              Filling the Element_Counter+2 : NZE  of the s_A and ij_A matrices with
              the off diagonal elements and their corresponding column indices
            */
            for (int index=0; index<13; index++)
              {
                /* Value of the Rest off-Diagonal Elements */
                s_A[Element_Counter +index] =
                  Special_Case.Value_Y[ Special_Case.Position_Y[index][1] ];

                /* Column index of the corresponding elements*/
                ij_A [Element_Counter + index] = Special_Case.Position_Y[index][0];
              }

            Element_Counter += 13;
           

          }

         if (j==ldy-1)
          {


            Special_Case.dy_special[0] = -1/
	    1./(2.*Delta_Y[k][j][i]
		*(Delta_Y[k][j-1][i]+Delta_Y[k][j][i]));
           
            Special_Case.dy_special[1] =1/
	    1./(2.*Delta_Y[k][j][i]
		*(Delta_Y[k][j-1][i]+Delta_Y[k][j][i]));
           
	    /*testing Purpuses */

		    // Special_Case.dy_special[0]=-1./(dy*dy);
		    // Special_Case.dy_special[1]=1./(dy*dy);



            for (int index=0; index<3; index++)
              {
                /*k-3:k-1*/
                Special_Case.Position_Y[index][0] = A(k_ref[index],j,i) +1;
                Special_Case.Position_Y[index][1] = index;
                Special_Case.Value_Y[index] = Coefficients_Z[index];

                /*j-1 outside loop*/

                /*i-3:i-1*/
                Special_Case.Position_Y[index+4][0] = A(k,j,i_ref[index]) +1;
                Special_Case.Position_Y[index+4][1] = index+4;
                Special_Case.Value_Y[index+4] = Coefficients_X[index];

                /*i+1:i+3*/
                Special_Case.Position_Y[index+7][0] = A(k,j,i_ref[index+3]) +1;
                Special_Case.Position_Y[index+7][1] = index+7;
                Special_Case.Value_Y[index+7] = Coefficients_X[2-index];

                /*k+1:k+3*/
                Special_Case.Position_Y[index+10][0] = A(k_ref[index+3],j,i) +1;
                Special_Case.Position_Y[index+10][1] = index+10;
                Special_Case.Value_Y[index+10] = Coefficients_Z[2-index];
              }


            /* j-1*/
            Special_Case.Position_Y[3][0] = A(k,j-1,i) +1;
            Special_Case.Position_Y[3][1] = 3;
            Special_Case.Value_Y[3] = Special_Case.dy_special[1];



            /* Diagonal Element */

            Special_Case.Diagonal = Coefficients_X[3] + Coefficients_Z[3]
              +Special_Case.dy_special[0];



            /* Calling the sorting function*/
            bubbleSort(13,Special_Case.Position_Y);



            /* Filling the first 1:Element_Counter elements of the s_A vector with the
               diagonal elements */
            s_A[A(k,j,i)+1] = Special_Case.Diagonal;

            /* Filling the Preconditioner vector with the diagonal Values of A*/
            Precond_A[A(k,j,i) + 1] =   Special_Case.Diagonal;


            /* Filling the first 1:Element_Counter elements of the ij_A vector with the
               position of the first off diagonal elements of each row*/
            ij_A [A(k,j,i)+1] = Element_Counter;


            /*
              Filling the Dimension_A+2 : NZE  of the s_A and ij_A matrices with
              the off diagonal elements and their corresponding column indices
            */



            for (int index=0; index<13; index++)
              {
                /* Value of the Rest off-Diagonal Elements */

                s_A[Element_Counter +index] =
                  Special_Case.Value_Y[ Special_Case.Position_Y[index][1] ];

                /* Column index of the corresponding elements*/
                ij_A [Element_Counter + index] = Special_Case.Position_Y[index][0];
              }

	    Element_Counter+=13;

	  }





      }
     
    }
  }



}


void bubbleSort( int size,int array[][2])
{
  int swapped;
  int i;
  for (i = 1; i < size; i++)
    {
      swapped = 0;    //this flag is to check if the array is already sorted
      int j;
      for(j = 0; j < size - i; j++)
        {
          if(array[j][0] > array[j+1][0])
            {
              int temp = array[j][0];
              int temp1 = array[j][1];
              array[j][0] = array[j+1][0];
              array[j][1] = array[j+1][1];
              array[j+1][0] = temp;
              array[j+1][1] = temp1;
              swapped = 1;
            }
        }
      if(!swapped){
        break; //if it is sorted then stop
      }
    }
}
