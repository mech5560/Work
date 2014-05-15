/*  Last Modified Time-stamp: <2014-05-01 18:15:27 mike_georgiou> */ 
//d Program that the problems constant's will be defined 
/* const identifier must be placed in order for the variables to be read only once */


// General contants

#ifndef Constants
#define Constants
const  double dt = 0.005;



// Non-Dimensional Numbers and quatities of my problem
const double Reynolds = 500.;
const double Pressure_Gradient = 12./Reynolds;
const double Prandtl  = 3.5;
const double Richardson  = 0.002;
#endif





#ifndef Boundary_Conditions
#define Boundary_Conditions

const double temperature_top = 1.;
const double temperature_bottom = 1.;

const double temperature_gradient_top = 0.;
const double temperature_gradient_bottom = 0.;

const double rho_top = 1.;
const double rho_bottom = 1.;

const double rho_gradient_top = 0.;
const double rho_gradient_bottom = 0.;

const double pressure_gradient_top = 0.;
const double pressure_gradient_bottom = 0.;

const double velocity_top=0.;
const double velocity_bottom=0.;

#endif
