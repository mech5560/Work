/*******************************************
 * Author: Michail Georgiou
 *  Last Modified: Time-stamp: <2014-04-30 11:36:29 mike_georgiou>
 *
 *
Output_Data.cpp -- This function will generate a .bin file containining
the data that are being passed to this function.

Input: The information that the user wants to send to the binary file
Output: The inputted information in binary format

*
* Written on Wednesday, 30 April 2014.
********************************************/

#include<fstream>
#include<sstream>

using namespace std;

void Output_Data(int ldx, int ldy, int ldz,
                 double*** dy, char* name)
{

	// info about stringstream
	// http://www.dreamincode.net/forums/topic/95826-stringstream-tutorial/
  stringstream output;
  
	string format = ".bin"; 

  output<<name<<format;

	//passing the stringstream content into another string
	const string tmp = output.str();
	// defining the filename of my output
	const char* filename = tmp.c_str();



  ofstream  outfile;
  outfile.open(filename,ios::out | ios::binary);


  for (int k=0; k<ldz; k++){
    for (int j=0; j<ldy; j++){
      for (int i=0; i<ldx; i++){

        outfile.write( (char*)&dy[k][j][i],sizeof(double));
 
     }
    }
  }



}
