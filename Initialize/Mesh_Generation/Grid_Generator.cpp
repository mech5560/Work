/*******************************************
 * Author: Michail Georgiou 
*  Last Modified: Time-stamp: <2014-04-30 11:12:40 mike_georgiou>   
*
*
Grid_Generator.cpp -- This program generates the coordinates for a
non-uniform mesh in the vertical direction. The results of this program will be
written in ".bin" file that at a later stage will be read by the main program 
in order to obtain the coordinates for the dy
*
* Written on Tuesday, 29 April 2014.
********************************************/
#include <iostream>
#include <cmath>
#include "Grid_Generator.h"

using namespace std;

int  main()
{

	cout<< "Please define the length in\nX (spanwise direction)\n";
	cout<<"Y (Vertical direction)\nZ (Streamwise direction)\n";

	//Reading the lengths given by the user
	double length_x,length_y,length_z;
	cin>>length_x; cin>>length_y;  cin>>length_z;

	cout<< "Please define the number of solution points in\nX (spanwise direction)\n";
	cout<<"Y (Vertical direction)\nZ (Streamwise direction)\n";


	//Reading the number of solution points in each direction.
	int points_x,points_y,points_z;
	cin>>points_x; cin>>points_y; 	cin>>points_z;


	//computing the dx and dz where they are constant for my case.
	double dx=length_x/points_x, dz=length_z/points_z;

	//Allocating memory for the dy matrix.
	double ***dy=Matrix_Alloc_New(points_z, points_y, points_x,
																0, 0,
																0, 0,
																0, 0);

	// cout<<"Please define if you want to generate uniform or non-uniform grid\n";
	// cout<<"For uniform press 0 and for non uniform 1\n";

	//Calling the function that will compute the hyperbolic mesh
	Hyperbolic_Mesh(points_x, points_y, points_z,
									dy);


	//Passing the data for the vertical directions to an output file
	char filename[]="dy_matrix";
	Output_Data(points_x, points_y, points_z, 
							dy,filename);



	//Releasing the allocated memory
	Free_Matrix_New(dy,
									0,0,0);
}
