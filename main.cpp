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
  length_x = 1.; length_y =1.; length_z=1.;
  ldx=64; ldy=64; ldz=64;

  //Calculating dx and dz
  dx= length_x/(ldx*1.0);
  dz= length_z/(ldz*1.0);

  dt= cfl*dx;

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


  Initial_Brown_2( Arr.velocity_x,  Arr.velocity_y,
  		   Arr.velocity_z,
  		   dx, Arr.dy,
  		   ldx,  ldy,  ldz);


  // Initial_Christos( Arr.velocity_x,  Arr.velocity_y,
  // 		    Arr.velocity_z,
  // 		    dx,
  // 		    ldx,  ldy,  ldz);

  BC_Velocities( Arr.velocity_x,
                 Arr.velocity_y,
                 Arr.velocity_z,
                 ldx, ldy, ldz,
                 left_x, right_x,
                 left_y, right_y,
                 left_z, right_z,
                 1,
                 velocity_x_top, velocity_x_bottom,
                 velocity_y_top, velocity_y_bottom,
                 velocity_z_top, velocity_z_bottom,
                 Arr.dy, dx, 0.);


  Print_2D_Matrix_Ghost(Arr.velocity_x, ldx,ldy,ldz,
			left_x, right_x,
			left_y, right_y,
			left_z, right_z,
			0,"Xinit");

  Print_2D_Matrix_Ghost(Arr.velocity_y, ldx,ldy,ldz,
			left_x, right_x,
			left_y, right_y,
			left_z, right_z,
			0,"Yinit");

  Print_2D_Matrix_Ghost(Arr.velocity_z, ldx,ldy,ldz,
			left_x, right_x,
			left_y, right_y,
			left_z, right_z,
			0,"Zinit");

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

  //////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

  //Only for the test case, this thing will be applied
  Initial_One(Arr.temperature_new,
              ldx, ldy, ldz,
              left_x,right_x,
              left_y,right_y,
              left_z,right_z);

  BC_Single(Arr.temperature_new,
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
            2,
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
            2,
            rho_gradient_top, rho_gradient_bottom,
            Arr.dy);



  //Only for the test case will thing will be applied
  //////////////////////////////////////////////////
  /////////////////////////////////////////////////
  Density_Calculator(Arr.rho_new,
                     Arr.temperature,
                     ldx, ldy, ldz);
  BC_Single(Arr.rho_new,
            ldx, ldy, ldz,
            left_x,right_x,
            left_y,right_y,
            left_z,right_z,
            2,
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

  //initializing the Arr.Residuals at the n-1 time
  //In order to do that I will use the Velosity Residual Functions.
  //but with the initial values velocities as an input

  double time_total =0.;
  Velocity_Residual_X( Arr.residual_x_old,
                       Arr.velocity_x, Arr.velocity_y,Arr.velocity_z,
                       Arr.flux_x,Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds,
                       Pressure_Gradient,
                       dx, Arr.dy,  dz,
                       time_total,
                       ldx,  ldy,  ldz);


  Velocity_Residual_Y( Arr.residual_y_old,
                       Arr.velocity_x,  Arr.velocity_y, Arr.velocity_z,
                       Arr.flux_x, Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds,
                       0.,
                       dx, Arr.dy,  dz,
                       time_total,
                       ldx,  ldy,  ldz);

  Velocity_Residual_Z( Arr.residual_z_old,
                       Arr.velocity_x,  Arr.velocity_y, Arr.velocity_z,
                       Arr.flux_x, Arr.flux_y,  Arr.flux_z,
                       Arr.temperature, Reynolds,
                       0.,
                       dx, Arr.dy, dz,
                       time_total,
                       ldx,  ldy, ldz);

  Print_2D_Matrix(Arr.residual_x_old, ldx,ldy,ldz,0,"rxold");
  Print_2D_Matrix(Arr.residual_y_old, ldx,ldy,ldz,0,"ryold");
  Print_2D_Matrix(Arr.residual_z_old, ldx,ldy,ldz,0,"rzold");



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
      time_total += dt;
      //////////////////////////////////////////////////////////////////////
      ///////////////////////////// Predictor Stage ////////////////////////
      //////////////////////////////////////////////////////////////////////

      // ////Calculating the Arr.temperature at the Predictor stage
      // Energy_Equation(Arr.temperature_new, Arr.temperature,
      //                 Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
      //                 Arr.rho, Reynolds, Prandtl,
      //                 dx, Arr.dy, dz, dt,
      //                 ldx, ldy, ldz);

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
      //          2,
      //           rho_gradient_top, rho_gradient_bottom,
      //           Arr.dy);


      /*Calculating the Intermediate Velocity at the Predictor stage*/
      Intermediate_Velocity_X(Arr.velocity_x_tilda,
                              Arr.residual_x, Arr.residual_x_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds,
                              Pressure_Gradient,
                              dx, Arr.dy, dz, dt,
                              time_total,
                              ldx,  ldy,  ldz);

      Intermediate_Velocity_Y(Arr.velocity_y_tilda,
                              Arr.residual_y, Arr.residual_y_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds,
                              0.,
                              dx, Arr.dy, dz, dt,
                              time_total,
                              ldx,  ldy,  ldz);

      Intermediate_Velocity_Z(Arr.velocity_z_tilda,
                              Arr.residual_z, Arr.residual_z_old,
                              Arr.velocity_x, Arr.velocity_y, Arr.velocity_z,
                              Arr.flux_x, Arr.flux_y, Arr.flux_z,
                              Arr.rho_new, Arr.rho, Arr.temperature_new,
                              Reynolds,
                              0.,
                              dx, Arr.dy, dz, dt,
                              time_total,
                              ldx,  ldy,  ldz);

      BC_Tilda(Arr.velocity_x_tilda,
               Arr.velocity_y_tilda,
               Arr.velocity_z_tilda,
               ldx, ldy, ldz,
               left_x, right_x,
               left_y, right_y,
               left_z, right_z);


      Print_2D_Matrix_Ghost(Arr.velocity_x_tilda, ldx,ldy,ldz,
			    left_x,right_x,
			    0,0,
			    0,0,
			    time_index,"tildax");
      Print_2D_Matrix_Ghost(Arr.velocity_y_tilda,
			    ldx,ldy,ldz,
			    0,0,
			    left_y,right_y,
			    0,0,
			    time_index,"tilday");


      Print_2D_Matrix_Ghost(Arr.velocity_z_tilda,
			    ldx,ldy,ldz,
			    0,0,
			    0,0,
			    left_z,right_z,
			    time_index,"tildaz");


      Print_2D_Matrix(Arr.residual_x, ldx,ldy,ldz,time_index,"resx");
      Print_2D_Matrix(Arr.residual_y, ldx,ldy,ldz,time_index,"resy");
      Print_2D_Matrix(Arr.residual_z, ldx,ldy,ldz,time_index,"resz");



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

      //Introducing as a  prediction the exact solution
      Print_1D_Matrix(rhs,
		  ldx,  ldy,  ldz,
		  time_index, "rhs");

      /*solving the Poisson Equation*/
      BCSG_Printing(s_a, ij_a, result, rhs, precond_a,
                    1e-15, 1000, dim_a, flag);
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
                2,
                pressure_gradient_top, pressure_gradient_bottom,
                Arr.dy);


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


      BC_Velocities( Arr.velocity_x_new,
                     Arr.velocity_y_new,
                     Arr.velocity_z_new,
                     ldx, ldy, ldz,
                     left_x, right_x,
                     left_y, right_y,
                     left_z, right_z,
                     1,
                     velocity_x_top, velocity_x_bottom,
                     velocity_y_top, velocity_y_bottom,
                     velocity_z_top, velocity_z_bottom,
                     Arr.dy, dx, time_total);


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


      /////////////////////////////////////////////////////////////
      ////////////////// end of the Predictor Stage ///////////////
      /////////////////////////////////////////////////////////////



      // /////////////////////////////////////////////////////////////
      // ////////////////// Corrector Stage //////////////////////////
      // /////////////////////////////////////////////////////////////


      // ////Calculating the Arr.temperature_Av which is necessary for the
      // //// calculation of the temperature
      // //// at the corrector stage
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
      //           2,
      //           rho_gradient_top, rho_gradient_bottom,
      //           Arr.dy);

      // /*Determining the RHS of the Poisson Equation*/
      // Right_Hand_Side_Poisson(rhs,Arr.velocity_x_tilda,
      //                         Arr.velocity_y_tilda,
      //                         Arr.velocity_z_tilda,
      //                         Arr.rho_new,Arr.rho, Arr.rho_old,
      //                         dx, Arr.dy, dz,
      //                         dt,
      //                         ldx,  ldy,  ldz);

      // /*Solving the Poisson Equation*/
      // BCSG_Printing(s_a, ij_a, result, rhs, precond_a,
      //               1e-15,1000,dim_a,flag);
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
      //           2,
      //           pressure_gradient_top, pressure_gradient_bottom,
      //           Arr.dy);

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
      // BC_Velocities( Arr.velocity_x_new,
      //                Arr.velocity_y_new,
      //                Arr.velocity_z_new,
      //                ldx, ldy, ldz,
      //                left_x, right_x,
      //                left_y, right_y,
      //                left_z, right_z,
      //                1,
      //                velocity_x_top, velocity_x_bottom,
      //                velocity_y_top, velocity_y_bottom,
      //                velocity_z_top, velocity_z_bottom,
      //                Arr.dy, dx, time_total);




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

      ////////////////////////////////////////////////////////////////
      ////////////////End of the Corrector Stage /////////////////////
      ///////////////////////////////////////////////////////////////



      if (time_index%1==0)
        {
          Print_2D_Matrix_Ghost(Arr.velocity_x_new,ldx,ldy,ldz,
				left_x,right_x,
				left_y,right_y,
				left_z,right_z,
				time_index,"X");

          Print_2D_Matrix_Ghost(Arr.velocity_y_new,ldx,ldy,ldz,
				left_x,right_x,
				left_y,right_y,
				left_z,right_z,
				time_index,"Y");

          Print_2D_Matrix_Ghost(Arr.velocity_z_new,ldx,ldy,ldz,
				left_x,right_x,
				left_y,right_y,
				left_z,right_z,
				time_index,"Z");

          Print_2D_Matrix_Ghost(Arr.pressure,ldx,ldy,ldz,
				left_x,right_x,
				left_y,right_y,
				left_z,right_z,
				time_index,"P");

          cout<<time_index<<endl;
          getchar();
        }

      Next_Step(&Arr);

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
