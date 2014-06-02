/*  Last Modified Time-stamp: <2014-05-27 14:27:59 mike_georgiou> */

#ifndef Functions
#define Functions

/* Allocating Functions*/
void Allocator( int ldx, int ldy, int ldz,
                int lx, int rx,
                int ly, int ry,
                int lz, int rz,
                Ar *Arr);

void DeAllocator( Ar *Arr,
                  int ldz, int ldy, int ldx,
                  int lx, int rx,
                  int ly, int ry,
                  int lz, int rz);

void Next_Step(Ar* Arr);



void Density_Calculator(double*** rho,
                        double*** T,
                        int ldx, int ldy, int ldz);



void BC_Corrector(Ar *Arr, int ldz, int ldy, int ldx,
                  int lz, int rz,
                  int ly, int ry,
                  int lx, int rx,
                  int index);

void BC_Predictor(Ar *Arr, int ldz, int ldy, int ldx,
                  int lz, int rz,
                  int ly, int ry,
                  int lx, int rx,
                  int index);




void BC_Single(double ***data,
               int ldx, int ldy, int ldz,
               int lx, int rx,
               int ly, int ry,
               int lz, int rz,
               int index,
               double bc_top, double bc_bottom,
               double *dy);



void BC_Tilda(double*** velocity_x_tilda, double*** velocity_y_tilda,
	      double*** velocity_z_tilda,
	      int ldx, int ldy, int ldz,
	      int lx, int rx,
	      int ly, int ry,
	      int lz, int rz);

void BC_Velocities(double*** velocity_x,
                   double*** velocity_y,
                   double*** velocity_z,
                   int ldx, int ldy, int ldz,
                   int lx, int rx,
                   int ly, int ry,
                   int lz, int rz,
                   int index,
                   double bc_x_top, double bc_x_bottom,
                   double bc_y_top, double bc_y_bottom,
                   double bc_z_top, double bc_z_bottom,
                   double *dy,double dx, double time);

/* Energy Equation Functions */
void Energy_Equation(double*** temperature_new, double*** temperature,
                     double*** velocity_x, double*** velocity_y, double***
                     velocity_z,
                     double*** rho, double Reynolds, double Prandtl,
                     double  dx, double* dy,double dz, double dt,
                     int ldx, int ldy,  int ldz);

void Energy_Equation_Corrector(double*** temperature_new,
                               double*** temperature_avg, double*** temperature,
                               double*** velocity_x, double*** velocity_y,
                               double*** velocity_z,
                               double*** rho, double Reynolds, double Prandtl,
                               double dx, double* dy, double dz, double dt,
                               int ldx, int ldy,  int ldz);

void Cubic_Mesh( double *dy,double *y,
		 double length_y,
                 int ldy,
                 int yl, int yr);

void Hyperbolic_Mesh( double *dy, double *y,
		      double length_y,
		      int ldy,
		      int yl, int yr);



void Initial_Zero( double ***Speed_X,
                   int ldz, int ldy, int ldx,
                   int lz, int rz,
                   int ly, int ry,
                   int lx, int rx);

void Initial_One( double ***Speed_X,
                  int ldz, int ldy, int ldx,
                  int lz, int rz,
                  int ly, int ry,
                  int lx, int rx);

void Initial_Reader(double*** velocity_x, double*** velocity_y,
		    double*** velocity_z,
		    int ldx, int ldy, int ldz);

void Initial_Christos(double*** velocity_x, double*** velocity_y,
		      double*** velocity_z,
		      double dx,
		      int ldx, int ldy, int ldz);

void Initial_Brown_2(double*** velocity_x, double*** velocity_y,
		     double*** velocity_z,
		     double dx, double* dy,
		     int ldx, int ldy, int ldz);

void Pertubation_Introducer(double*** velocity_x, double***velocity_y, 
			    double*** velocity_z,
			    double length_x, double length_y, 
			    double length_z,
			    double dx, double *dy, double dz,
			    int ldx, int ldy, int ldz);



void Initial_Conditions_Turbulence(double*** velocity_x,
                                   double*** velocity_y,
                                   double*** velocity_z,
                                   double Reynolds, double dt,
                                   double dx, double* dy, double dz,
                                   int ldx, int ldy, int ldz);


void Print_3D(double ***A,
	      int ldz, int ldy, int ldx,
	      int time_index, char *mike);

void Print_2D_Curve(double*** data, double* dy,
                   int ldx, int ldy, int ldz,
                   int time_index, char* name);


void Print_2D_Matrix(double ***A,
		     int ldz, int ldy, int ldx,
		     int time_index, char*);

void Print_2D_Matrix_Ghost(double ***A,
			   int ldz, int ldy, int ldx,
			   int, int,
			   int, int,
			   int, int,
			   int time_index, char*);


void Print_1D_Matrix(double *A,
		     int ldz, int ldy, int ldx,
		     int time_index, char*);


void Print_3D_Binary(double ***A,
                     int ldz, int ldy, int ldx,
                     int time_index, char *mike);


void Velocity_Residual_X( double*** residual_x, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y, double*** flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
			  double time,
                          int ldx, int ldy, int ldz);

void Velocity_Residual_Y( double*** residual_y,double*** velocity_x,
                          double*** velocity_y,double*** velocity_z,
                          double*** flux_x,double*** flux_y, double***flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
			  double time,
                          int ldx, int ldy, int ldz);


void Velocity_Residual_Z( double*** residual_z, double*** velocity_x,
                          double*** velocity_y, double*** velocity_z,
                          double*** flux_x, double*** flux_y,
			  double*** flux_z,
                          double*** temperature, double Reynolds,
			  double source,
                          double dx, double* dy, double dz,
			  double time,
                          int ldx, int ldy, int ldz);




void Intermediate_Velocity_X(double*** velocity_x_tilda,
                             double*** residual_x, double*** residual_x_old,
                             double*** velocity_x, double*** velocity_y,
                             double*** velocity_z,
                             double*** flux_x, double***flux_y, double***flux_z,
                             double*** rho_new, double*** rho,
                             double***  temperature, double Reynolds,
			     double source,
                             double dz, double* dy,  double dx,
			     double dt, double time_total,
                             int ldx, int ldy, int ldz);

void Intermediate_Velocity_Y(double*** velocity_y_tilda,
                             double*** residual_y, double*** residual_y_old,
                             double*** velocity_x, double*** velocity_y,
                             double*** velocity_z,
                             double*** flux_x, double***flux_y, double***flux_z,
                             double*** rho_new, double*** rho,
                             double*** temperature, double  Reynolds,
			     double source,
                             double dx, double* dy,  double dz,
			     double dt, double time_total,
                             int ldx, int ldy, int ldz);

void Intermediate_Velocity_Z(double*** velocity_z_tilda,
                             double*** residual_z, double*** residual_z_old,
                             double*** velocity_x, double*** velocity_y,
                             double*** velocity_z,
                             double*** flux_x, double***flux_y, double***flux_z,
                             double*** rho_new, double*** rho,
                             double*** temperature, double Reynolds, double source,
                             double dx, double* dy,  double dz, 
			     double dt, double time_total,
                             int ldx, int ldy, int ldz);


/* Poisson Solver*/
void BCSG(double *sL, int *ijL, double *X, double *D, double *Pre,
          double error, int L, int N, int& flag);

void BCSG_Printing(double *sL, int *ijL, double *X, double *D, double *Pre,
                   double error, int L, int N, int& flag);


void Vector_Constructor(double *s_A, double *Precond_A, int *ij_A,
                        double Coefficients_Z[4], double Coefficients_X[4], double *dy,
                        int ldz, int ldy, int ldx, int Nze);


void Right_Hand_Side_Poisson(double* rhs, double*** velocity_x,
                             double*** velocity_y, double*** velocity_z,
                             double*** rho_new, double*** rho,double*** rho_old,
                             double dx, double* dy, double dz,
                             double dt,
                             int ldx, int ldy, int ldz);

void RHS_Poisson( double *RHS, double ***Speed_Int_X,
                  double ***Speed_Int_Y, double ***Speed_Int_Z,
                  double ***Rho_Star, double ***Rho_N, double ***Rho_Nm1, double *dy,
                  double dz, double dx,int ldz, int ldy, int ldx);



/* Velocity Update*/
void Velocity_Update_X(double*** velocity_x, double*** velocity_x_tilda,
                       double*** rho, double*** pressure,
                       double dx, double dt,
                       int ldx, int ldy, int ldz);


void Velocity_Update_Y(double*** velocity_y, double*** velocity_y_tilda,
                       double*** rho, double*** pressure,
                       double* dy, double dt,
                       int ldx, int ldy, int ldz);

void Velocity_Update_Z(double*** velocity_z, double*** velocity_z_tilda,
                       double*** rho, double*** pressure,
                       double dz, double dt,
                       int ldx, int ldy, int ldz);




/* Flux Evaluation */
void Flux_Evaluation_X(double*** flux_x, double*** velocity_x,
                       double*** rho, double*** pressure,
                       double dx, double dt,
                       int ldx, int ldy, int ldz);

void Flux_Evaluation_Y(double*** flux_y, double*** velocity_y,
                       double*** rho, double*** pressure,
                       double* dy, double dt,
                       int ldx, int ldy, int ldz);

void Flux_Evaluation_Z(double*** flux_z, double*** velocity_z,
                       double*** rho, double*** pressure,
                       double dz, double dt,
                       int ldx, int ldy, int ldz);


void Intermediate_Velocity_Z_Press(double*** velocity_z_tilda,
                                   double*** residual_z,
                                   double*** residual_z_old,
                                   double*** velocity_x,
                                   double*** velocity_y,
                                   double*** velocity_z,
                                   double*** flux_x, double***flux_y,
                                   double***flux_z,
                                   double*** rho_new, double*** rho,
                                   double*** pressure,
                                   double*** temperature,
                                   double Reynolds, double source_term,
                                   double dx, double* dy,  double dz,
                                   double dt, double time_total,
                                   int ldx, int ldy, int ldz);

void Intermediate_Velocity_Y_Press(double*** velocity_y_tilda,
				   double*** residual_y,
				   double*** residual_y_old,
				   double*** velocity_x, 
				   double*** velocity_y,
				   double*** velocity_z,
				   double*** flux_x, double***flux_y, 
				   double***flux_z,
				   double*** rho_new, double*** rho,
				   double*** temperature,
				   double*** pressure,
				   double Reynolds,double source,
				   double dx, double* dy,  double dz,
				   double dt, double time_total,
				   int ldx, int ldy, int ldz);


void Intermediate_Velocity_X_Press(double*** velocity_x_tilda,
				   double*** residual_x,
				   double*** residual_x_old,
				   double*** velocity_x, 
				   double*** velocity_y,
				   double*** velocity_z,
				   double*** flux_x, double***flux_y,
				   double***flux_z,
				   double*** rho_new, double*** rho,
				   double*** pressure,
				   double*** temperature,
				   double Reynolds,double source,
				   double dx, double* dy,  double dz,
				   double dt, double time_total,
				   int ldx, int ldy, int ldz);



void Continuity_Printer_Y(double*** velocity_y, double* dy,
                          int ldx, int ldy, int ldz,
                          int time_index, char *mike);


#endif
