/*******************************************
 *   Author: Michail Georgiou
 * main.cpp -- This project solves the 3D navier stokes equations for a
non-isothermal problem. Fourth order accurate conservative schemes are used in
the spanwise(x) and streamwise(z) directions. These schemes are being adopted by
the Lessani-Papalexandris paper. On the vertical direction(y) a non-uniform grid
is used.
*
* Written on Wednesday, 30 April 2014.
********************************************/

//Definition of the libraries that this program will use
#include "main.h"
#include "main-inl.h"

#include "./Header_Files/Functions.h"
#include "./Header_Files/Data.h"

int main (int argc, char *argv[])
{

  //Defining the number of solution points in each direction
  length_x = 4.; length_y =2.; length_z=2.;

  ldx=62; ldy=31; ldz=31;

  //Calculating dx and dz
  dx= length_x/(ldx*1.0);
  dz= length_z/(ldz*1.0);

  dt= cfl*dz;

  //Defining the Non-Zero Elements for the Poisson solver
  nze = 15*ldz*(ldy-2)*ldx + 14*ldz*ldx*2;

  //Dimension of the A matrix of the Poisson Equation
  dim_a=ldz*ldy*ldx;

  //Defining the extra points on each direction depending on the BC
  left_x=3; right_x=3;
  left_y=1; right_y=1;
  left_z=3; right_z=3;

  Allocator(ldx,ldy,ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            &Arr);


  //Vectors for the Poisson Solver
  s_a = new double [nze  +2];
  ij_a = new int [nze  +2];
  rhs = new double [dim_a +1];
  precond_a = new double [dim_a +1];
  result = new double [dim_a +1];

  //Calculation of the Dy's
  Cubic_Mesh(Arr.dy, Arr.y,
	     length_y,
	     ldy,
	     left_y,right_y);


  // Defining the Second derivative coefficients  that will be use by the
  // Vector_Constructor Function
  Coefficients_X[0] = 1./(576.*dx*dx);
  Coefficients_X[1] = -9./(96.*dx*dx);
  Coefficients_X[2] = ( 81./(64.*dx*dx) + 9./(96.*dx*dx) );
  Coefficients_X[3] = (- 81./(32.*dx*dx) - 1./(288.*dx*dx) );

  Coefficients_Z[0] = 1./(576.*dz*dz);
  Coefficients_Z[1] = -9./(96.*dz*dz);
  Coefficients_Z[2] = ( 81./(64.*dz*dz) + 9./(96.*dz*dz) );
  Coefficients_Z[3] = (- 81./(32.*dz*dz) - 1./(288.*dz*dz) );


  /////////////////////////////////////////////////////////////
  //////////////////Initial Conditions/////////////////////////
  /////////////////////////////////////////////////////////////

  //Assigning zero Initial conditions for the velocities and Arr.Temperature

  ////Initializing the Velocities
  Initial_Reader( Arr.velocity_x,  Arr.velocity_y,
  		  Arr.velocity_z,
  		  ldx,  ldy,  ldz);

  Pertubation_Introducer( Arr.velocity_x,  Arr.velocity_y,
  			  Arr.velocity_z,
  			  length_x,  length_y,
  			  length_z,
  			  dx,  Arr.y,  dz,
  			  ldx,  ldy,  ldz);

  // Print_2D_Curve(Arr.velocity_x, Arr.dy,ldz,ldy,ldx,0,"X0");
  // Print_2D_Curve(Arr.velocity_y, Arr.dy,ldz,ldy,ldx,0,"Y0");

  BC_Single(Arr.velocity_x,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            1,
            velocity_top, velocity_bottom,
            Arr.dy);

  BC_Single(Arr.velocity_y,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            1,
            velocity_top, velocity_bottom,
            Arr.dy);

  BC_Single(Arr.velocity_z,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            1,
            velocity_top, velocity_bottom,
            Arr.dy);

  //Initializing the Arr.temperature
  Initial_One(Arr.temperature,
              ldx, ldy, ldz,
              left_x,right_x,
              left_y,right_y,
              left_z,right_z);


  BC_Single(Arr.temperature,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            1,
            temperature_top, temperature_bottom,
            Arr.dy);


  //Initializing  the Density for each Point
  Density_Calculator(Arr.rho_old,
                     Arr.temperature,
                     ldx, ldy, ldz);

  BC_Single(Arr.rho_old,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            0,
            rho_gradient_top, rho_gradient_bottom,
            Arr.dy);


  Density_Calculator(Arr.rho,
                     Arr.temperature,
                     ldx, ldy, ldz);

  BC_Single(Arr.rho,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
	    left_z,right_z,
            0,
            rho_gradient_top, rho_gradient_bottom,
            Arr.dy);


  //Initializing the Arr.Fluxes.
  //In order to do that I will use the Arr.Flux_Evaluation function.
  //But with the initial values as an input

  Flux_Evaluation_X(Arr.flux_x, Arr.velocity_x,
                    Arr.rho_old, Arr.pressure,
                    dx, dt,
                    ldx+1,  ldy,  ldz);

  Flux_Evaluation_Y(Arr.flux_y, Arr.velocity_y,
                    Arr.rho_old, Arr.pressure,
                    Arr.dy, dt,
                    ldx,  ldy+1, ldz);

  Flux_Evaluation_Z(Arr.flux_z, Arr.velocity_z,
                    Arr.rho_old,  Arr.pressure,
                    dz, dt,
                    ldx,  ldy,  ldz+1);

  //Assigning the BC for the Fluxes and the Intermediate Velocities
  BC_Flux(Arr.flux_x, Arr.flux_y, Arr.flux_z,
          Arr.velocity_x_tilda, Arr.velocity_y_tilda, Arr.velocity_z_tilda,
          ldx, ldy, ldz,
          left_x, right_x,
          left_y, right_y,
          left_z, right_z);


  //Initializing the Arr.Residuals at the n-1 time
  //In order to do that I will use the Velosity Residual Functions.
  //but with the initial values velocities as an input
  Velocity_Residual_X( Arr.residual_x_old,
                       Arr.velocity_x, Arr.velocity_y,Arr.velocity_z,
                       Arr.flux_x,Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds,
		       Pressure_Gradient,
                       dx, Arr.dy,  dz,
                       ldx,  ldy,  ldz);


  Velocity_Residual_Y( Arr.residual_y_old,
                       Arr.velocity_x,  Arr.velocity_y, Arr.velocity_z,
                       Arr.flux_x, Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds,
		       0.,
                       dx, Arr.dy,  dz,
                       ldx,  ldy,  ldz);

  Velocity_Residual_Z( Arr.residual_z_old,
                       Arr.velocity_x,  Arr.velocity_y, Arr.velocity_z,
                       Arr.flux_x, Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds, 
		       0.,
                       dx, Arr.dy, dz,
                       ldx,  ldy, ldz);

  // Before Entering the time integration loop
  // i will define the constant vectors for the poisson solver
  Vector_Constructor(s_a, precond_a, ij_a,
                     Coefficients_Z,  Coefficients_X, Arr.dy,
                     ldx,  ldy,  ldz,  nze);

  //////////////////////////////////////////////////////////////////////
  //////////////////////Time Integration Loop///////////////////////////
  //////////////////////////////////////////////////////////////////////
  for (int time_index=0; time_index<4e5; time_index++ )
    {

      //////////////////////////////////////////////////////////////////////
      ///////////////////////////// Predictor Stage ////////////////////////
      //////////////////////////////////////////////////////////////////////

      //Calculating the Arr.temperature at the Predictor stage
      Energy_Equation(Arr.temperature_new, Arr.temperature,
                      Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                      Arr.rho, Reynolds, Prandtl,
                      dx, Arr.dy, dz, dt,
                      ldx, ldy, ldz);

      BC_Single(Arr.temperature_new,
                ldx, ldy, ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                1,
                temperature_top, temperature_bottom,
                Arr.dy);


      /*Computing Arr.rho_Star*/
      Density_Calculator(Arr.rho_new, Arr.temperature_new, ldx, ldy, ldz);
      BC_Single(Arr.rho_new,
                ldx,ldy,ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                0,
                rho_gradient_top, rho_gradient_bottom,
                Arr.dy);


      /*Calculating the Intermediate Velocity at the Predictor stage*/
      Intermediate_Velocity_X(Arr.velocity_x_tilda,
                              Arr.residual_x, Arr.residual_x_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds,
			      Pressure_Gradient,
                              dx, Arr.dy, dz, dt,
                              ldx,  ldy,  ldz);

      Intermediate_Velocity_Y(Arr.velocity_y_tilda,
                              Arr.residual_y, Arr.residual_y_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds, 0.,
                              dx, Arr.dy, dz, dt,
                              ldx,  ldy,  ldz);

      Intermediate_Velocity_Z(Arr.velocity_z_tilda,
                              Arr.residual_z, Arr.residual_z_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds, 
			      0.,
                              dx, Arr.dy, dz, dt,
                              ldx,  ldy,  ldz);


      BC_Flux(Arr.flux_x, Arr.flux_y, Arr.flux_z,
              Arr.velocity_x_tilda, Arr.velocity_y_tilda, Arr.velocity_z_tilda,
              ldx, ldy, ldz,
              left_x, right_x,
              left_y, right_y,
              left_z, right_z);


      /*Solving the Poisson Equation
        To do that, I will use the bcgs-solver of Christos.
        The s_A ij_A and Precond_A vectors will be constructed once and
        they will be defined outside the iterating part of the code*/


      //Determining the RHS of the Poisson Equation
      Right_Hand_Side_Poisson(rhs,Arr.velocity_x_tilda,
                              Arr.velocity_y_tilda, Arr.velocity_z_tilda,
                              Arr.rho_new,Arr.rho, Arr.rho_old,
                              dx, Arr.dy, dz,
                              dt,
                              ldx, ldy,  ldz);

      /*solving the Poisson Equation*/
      BCSG(s_a, ij_a, result, rhs, precond_a,
		    1e-10, 3000, dim_a, flag);
      if(flag==1)
        {
          cout<<time_index<<endl;
          break;
        }

      /*Passing the Poisson solution to the Arr.pressure Array*/
      for (int k=0; k<ldz; k++){
        for (int j=0; j<ldy; j++){
          for (int i=0; i<ldx; i++){

            Arr.pressure[k][j][i] = result[A(i,ldx,j,ldy,k,ldz) +1];
          }
        }
      }
      BC_Single(Arr.pressure,
                ldx,ldy,ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                0,
                pressure_gradient_top, pressure_gradient_bottom,
                Arr.dy);

      //Print_2D_Matrix(Arr.pressure,ldz,ldy,ldx,time_index,"prep");
      
      /*computing the Updated Velocity*/
      Velocity_Update_X(Arr.velocity_x_new, Arr.velocity_x_tilda,
                        Arr.rho_new, Arr.pressure,
                        dx, dt,
                        ldx,  ldy,  ldz);

      Velocity_Update_Y(Arr.velocity_y_new,  Arr.velocity_y_tilda,
                        Arr.rho_new, Arr.pressure,
                        Arr.dy, dt,
                        ldx,  ldy,  ldz);


      Velocity_Update_Z(Arr.velocity_z_new,  Arr.velocity_z_tilda,
                        Arr.rho_new, Arr.pressure,
                        dz, dt,
                        ldx,  ldy,  ldz);

      BC_Single(Arr.velocity_x_new,
                ldx, ldy, ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                1,
                velocity_top, velocity_bottom,
                Arr.dy);

      BC_Single(Arr.velocity_y_new,
                ldx, ldy, ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                1,
                velocity_top, velocity_bottom,
                Arr.dy);

      BC_Single(Arr.velocity_z_new,
                ldx, ldy, ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z,
                1,
                velocity_top, velocity_bottom,
                Arr.dy);

      /*Updating the Auxiliary Arr.Fluxes in order to proceed at
        the Corrector Stage
      */
      Flux_Evaluation_X(Arr.flux_x, Arr.velocity_x_tilda,
                        Arr.rho_new, Arr.pressure,
                        dx, dt,
                        ldx+1,  ldy,  ldz);


      Flux_Evaluation_Y(Arr.flux_y, Arr.velocity_y_tilda,
                        Arr.rho_new, Arr.pressure,
                        Arr.dy, dt,
                        ldx,  ldy+1,  ldz);


      Flux_Evaluation_Z(Arr.flux_z, Arr.velocity_z_tilda,
                        Arr.rho_new,  Arr.pressure,
                        dx, dt,
                        ldx,  ldy,  ldz+1);



      /*Applying the Boundary conditions to the computed fluxes*/
      BC_Flux(Arr.flux_x, Arr.flux_y, Arr.flux_z,
              Arr.velocity_x_tilda, Arr.velocity_y_tilda, Arr.velocity_z_tilda,
              ldx, ldy, ldz,
              left_x, right_x,
              left_y, right_y,
              left_z, right_z);


      /////////////////////////////////////////////////////////////
      ////////////////// end of the Predictor Stage ///////////////
      /////////////////////////////////////////////////////////////



      // /////////////////////////////////////////////////////////////
      // ////////////////// Corrector Stage //////////////////////////
      // /////////////////////////////////////////////////////////////


      // //Calculating the Arr.temperature_Av which is necessary for the
      // // calculation of the temperature 
      // // at the corrector stage
      // for (int k=0; k<ldz; k++){
      //   for (int j =0; j<ldy; j++){
      //     for (int i=0; i<ldx; i++){

      //       Arr.temperature_avg[k][j][i]=0.5*
      //         (Arr.temperature[k][j][i] + Arr.temperature_new[k][j][i]);

      //     }
      //   }
      // }
      // BC_Single(Arr.temperature_avg,
      //           ldx,ldy,ldz,
      //           left_z,right_z,
      //           left_y,right_y,
      //           left_x,right_x,1,
      //           temperature_top, temperature_bottom, Arr.dy);


      // /*Calculating the Arr.temperature at the Predictor stage*/
      // Energy_Equation_Corrector(Arr.temperature_new,Arr.temperature_avg,
      //                           Arr.temperature,
      //                           Arr.velocity_z_new, Arr.velocity_y_new,
      //                           Arr.velocity_x_new,
      //                           Arr.rho_new, Reynolds, Prandtl,
      //                           dx, Arr.dy, dz, dt,
      //                           ldx,  ldy,   ldz);

      // BC_Single(Arr.temperature_new,
      //           ldx, ldy, ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           1,
      //           temperature_top, temperature_bottom,
      //           Arr.dy);


      // /*Computing Arr.rho_Star*/
      // Density_Calculator(Arr.rho_new, Arr.temperature_new, ldx, ldy, ldz);
      // BC_Single(Arr.rho_new,
      //           ldx,ldy,ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           0,
      //           rho_gradient_top, rho_gradient_bottom,
      //           Arr.dy);

      // /*Determining the RHS of the Poisson Equation*/
      // Right_Hand_Side_Poisson(rhs,Arr.velocity_x_tilda,
      //                         Arr.velocity_y_tilda, Arr.velocity_z_tilda,
      //                         Arr.rho_new,Arr.rho, Arr.rho_old,
      //                         dx, Arr.dy, dz,
      //                         dt,
      //                         ldx,  ldy,  ldz);

      // /*Solving the Poisson Equation*/
      // BCSG(s_a, ij_a, result, rhs, precond_a,
      // 		    1e-10,3000,dim_a,flag);
      // if(flag==1)
      //   {
      //     cout<<time_index<<endl;
      //     break;
      //   }

      // /*Passing the Poisson solution to the Arr.pressure Array*/
      // for (int k=0; k<ldz; k++){
      //   for (int j=0; j<ldy; j++){
      //     for (int i=0; i<ldx; i++){

      //       Arr.pressure[k][j][i] = result[A(i,ldx,j,ldy,k,ldz) +1];

      //     }
      //   }
      // }

      // BC_Single(Arr.pressure,
      //           ldx,ldy,ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           0,
      //           pressure_gradient_top, pressure_gradient_bottom,
      //           Arr.dy);

      // Print_2D_Matrix(Arr.pressure,ldz,ldy,ldx,time_index,"prec");


      // /*computing the Updated Velocity*/
      // Velocity_Update_X(Arr.velocity_x_new, Arr.velocity_x_tilda,
      //                   Arr.rho_new, Arr.pressure,
      //                   dx, dt,
      //                   ldx,  ldy,  ldz);

      // Velocity_Update_Y(Arr.velocity_y_new,  Arr.velocity_y_tilda,
      //                   Arr.rho_new, Arr.pressure,
      //                   Arr.dy, dt,
      //                   ldx,  ldy,  ldz);


      // Velocity_Update_Z(Arr.velocity_z_new,  Arr.velocity_z_tilda,
      //                   Arr.rho_new, Arr.pressure,
      //                   dz, dt,
      //                   ldx,  ldy,  ldz);
      // /*Implementing the Velocity  boundary conditions for the next step*/

      // BC_Single(Arr.velocity_x,
      //           ldx, ldy, ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           1,
      //           velocity_top, velocity_bottom,
      //           Arr.dy);

      // BC_Single(Arr.velocity_y,
      //           ldx, ldy, ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           1,
      //           velocity_top, velocity_bottom,
      //           Arr.dy);

      // BC_Single(Arr.velocity_z,
      //           ldx, ldy, ldz,
      //           left_x,right_x,
      //           left_y,right_y,
      //           left_z,right_z,
      //           1,
      //           velocity_top, velocity_bottom,
      //           Arr.dy);

      // /*Updating the Auxiliary Arr.Fluxes in order to proceed at
      //   the next timestep*/
      // Flux_Evaluation_X(Arr.flux_x, Arr.velocity_x_tilda,
      //                   Arr.rho_new, Arr.pressure,
      //                   dx,dt,
      //                   ldx+1,  ldy,  ldz);


      // Flux_Evaluation_Y(Arr.flux_y, Arr.velocity_y_tilda,
      //                   Arr.rho_new, Arr.pressure,
      //                   Arr.dy, dt,
      //                   ldx,  ldy+1,  ldz);


      // Flux_Evaluation_Z(Arr.flux_z, Arr.velocity_z_tilda,
      //                   Arr.rho_new,  Arr.pressure,
      //                   dx, dt,
      //                   ldx,  ldy,  ldz+1);



      // BC_Flux(Arr.flux_x, Arr.flux_y, Arr.flux_z,
      //         Arr.velocity_x_tilda, Arr.velocity_y_tilda, Arr.velocity_z_tilda,
      //         ldx, ldy, ldz,
      //         left_x, right_x,
      //         left_y, right_y,
      //         left_z, right_z);

      ////////////////////////////////////////////////////////////////
      ////////////////End of the Corrector Stage /////////////////////
      ///////////////////////////////////////////////////////////////

      // Enforcing the Z component of the velocity to be zero
      Initial_Zero(Arr.velocity_z_new,
		  ldx, ldy, ldz,
		  left_x,right_x,
		  left_y,right_y,
		  left_z,right_z);




      if (time_index%100==0)
        {
          Print_2D_Matrix(Arr.velocity_x_new,ldz,ldy,ldx,time_index,"X"); 
          Print_2D_Matrix(Arr.velocity_y_new,ldz,ldy,ldx,time_index,"Y");
          Print_2D_Matrix(Arr.velocity_z_new,ldz,ldy,ldx,time_index,"Z");
          Print_2D_Matrix(Arr.pressure,ldz,ldy,ldx,time_index,"P");

          Print_2D_Curve(Arr.velocity_x_new, Arr.dy,ldz,ldy,ldx,time_index,"2d");
          cout<<time_index<<endl;

        }


      /*Passing the data to proceed  to the next time step*/
      Temp=Arr.velocity_x;
      Arr.velocity_x= Arr.velocity_x_new;
      Arr.velocity_x_new = Temp;

      Temp=Arr.velocity_y;
      Arr.velocity_y= Arr.velocity_y_new;
      Arr.velocity_y_new = Temp;

      Temp=Arr.velocity_z;
      Arr.velocity_z= Arr.velocity_z_new;
      Arr.velocity_z_new = Temp;

      Temp=Arr.temperature;
      Arr.temperature = Arr.temperature_new;
      Arr.temperature_new = Temp;

      Temp = Arr.rho_old;
      Arr.rho_old = Arr.rho;
      Arr.rho = Temp;

      Temp = Arr.rho;
      Arr.rho = Arr.rho_new;
      Arr.rho_new = Temp;

      Temp = Arr.residual_x_old;
      Arr.residual_x_old = Arr.residual_x;
      Arr.residual_x = Temp;

      Temp = Arr.residual_y_old;
      Arr.residual_y_old = Arr.residual_y;
      Arr.residual_y = Temp;

      Temp = Arr.residual_z_old;
      Arr.residual_z_old = Arr.residual_z;
      Arr.residual_z = Temp;
    }

  /*Releasing the allocated memory*/
  DeAllocator(  &Arr,
                ldx,  ldy,  ldz,
                left_x,right_x,
                left_y,right_y,
                left_z,right_z);

  delete[] s_a;
  delete[] ij_a;
  delete[] precond_a;
  delete[] rhs;
  delete[] result;
}
