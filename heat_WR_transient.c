// Heat equation in 2d with random numbers for steady state
//  gcc -o heat_WR_transient.exe heat_WR_transient.c  ranlib.c rnglib.c -fopenmp
// ./heat_WR_transient.exe
// Main parameters
#define N_WALKERS 0 // number of parallel walkers 0 for maximum hardware number of threads
#define MAX_VARIANCE_TO_MEAN_RATIO -1.0 // maximum variance to mean ratio to stop the simulation (-1 for no control by variance to ratio)
#define DEEP_OF_DOMAIN 10.0 // deep of domain in meters
#define WIDTH_OF_DOMAIN 10.0 // width of domain in meters
#define NUMBER_OF_POINTS_IN_X 10000 // 7500
#define NUMBER_OF_POINTS_IN_Y 10000 //7500
#define MAX_NUMBER_OF_LOOPS 10000 // maximum number of loops
#define DEEP_CABLE 1.04 // deep of cable in meters
#define SEP_CABLES 0.4  // separation between cables in meters between three phase circuits
#define TYPE_CABLE_ARRANGEMENT 2 // 1 for horizontal arrangement, 2 for trifoil arrangement
#define BACKFILL_HEIGHT_LAYER 0.4 //1.04 // 1.04 // height of backfill soil layer in m
#define VERBOSE 0 // 1 to print detailed results of total results in the screen
#define VERBOSE2 0 // 1 to print detailed results of total results in the screen
#define VERBOSE3 0 // 1 to print detailed results of probabilities in the screen
#define VERBOSE4 0 // 1 to print detailed results of heat sources in the screen
#define PRINT_ADVANCE 1 // 1 to print advance of the simulation in the screen
#define RANDOM_LIB 1 // 1 to use native random library, 0 to use ranlib library
#define TMAX 36000.0 // 360000.0 // maximum temperature of the domain 
#define Factor_dt 1.0 // factor to calculate the time step

// libraries
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
# include "ranlib.h"
# include "rnglib.h"

// Global Variables

// Main parameters
double L = DEEP_OF_DOMAIN; // length of domain
double W = WIDTH_OF_DOMAIN; // width of domain
int nx= NUMBER_OF_POINTS_IN_X; // number of points in x direction
int ny= NUMBER_OF_POINTS_IN_Y; // number of points in y direction

double hmin; // spatial step
int n_steps; // number of steps in each walk
double d_t_0; // time step

// Aluminum 630 mm2 data for calculations
double cable_area= 0.000630;// cross area of cable 630 mm2 in squered meters
double cable_diam = 0.0615; // diameter of cable 630 mm2 - fuente: http://es.lemeicable.com/product/108.html
double cond_diam =  0.0301; // diameter of conductor 630 mm2 - fuente: http://es.lemeicable.com/product/108.html 
double res_cond =   0.000064; // ohm/m  - electric resistivity of cable 630 mm2 0.0469 ohm/km at 20°C DC - fuente: http://es.lemeicable.com/product/108.html
double rho_cond = 2712; //  cobre->8960.0; // Aluminium density  in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_cond = 897.0; // specific heat of copper cable in J/(kg K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double k_cond = 237.0; // thermal conductivity of aluminum cable in W/(m.K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html

// datos XLPE
double rho_xlpe = 920.0; // density of XLPE in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_xlpe = 2174.0; // specific heat of XLPE in J/(kg K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
double k_xlpe = 0.28; // thermal conductivity of XLPE in W/(m.K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
double F_xlpe = 0.0; // heat source in cable W/m3

// Native Soil thermal data for calculations
double rho_soil1 = 2000.0; // density of soil_1 in kg/m3 - fuente: https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil1 = 800.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil1 = (1.0/1.2) ;// 0.2; // thermal conductivity of soil 1 in W/(m.K) (resistivity = 5)- fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5
double F_soil1 = 0.0; // heat source in soil 1 W/m3

// Backfill Soil thermal data for calculations
double rho_soil2 = 1500.0; // density of soil_1 in kg/m3 - fuente:  https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil2 = 900.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil2 = (1.0/1.2); // 0.3333; // thermal conductivity of soil 1 in W/(m.K) (resistivity = 3) - fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5 
double soil2_layer = BACKFILL_HEIGHT_LAYER; // height of backfill soil layer in m
double F_soil2 = 0.0; // heat source in soil 1 W/m3

// max, min values
double rho_min = 1500; // minimum density
double Ce_min = 900; // minimum specific heat
double k_max = 237; // maximum thermal diffusivity

//Location of conductors (x=0 is the left side of the domain, y=0 is the top of the domain)
double deep_cable = DEEP_CABLE; // depth of cable in m, from the surface to the cable touch the soil
double separation_circuits = SEP_CABLES; // separation between cable circuits in m

//Circuit 1 - three phase arrangement in trifoil - Central circuit
double x_cable11; // x coordinate of cable 1 , central and upper cable
double y_cable11; // y coordinate of cable 1
double x_cable12; // x coordinate of cable 2, left and down cable
double y_cable12; // y coordinate of cable 2
double x_cable13; // x coordinate of cable 3, right and down cable
double y_cable13; // y coordinate of cable 3
double F_c1; // heat source in cable W/m3

//Circuit 2 - three phase arrangement in trifoil - Left circuit
double x_cable21; // x coordinate of cable 1 , central and upper cable
double y_cable21; // y coordinate of cable 1
double x_cable22; // x coordinate of cable 2, left and down cable
double y_cable22; // y coordinate of cable 2
double x_cable23; // x coordinate of cable 3, right and down cable
double y_cable23; // y coordinate of cable 3
double F_c2; // heat source in cable W/m3

// //Circuit 3 - three phase arrangement in trifoil - Right circuit
// double x_cable31; // x coordinate of cable 1 , central and upper cable
// double y_cable31; // y coordinate of cable 1
// double x_cable32; // x coordinate of cable 2, left and down cable
// double y_cable32; // y coordinate of cable 2
// double x_cable33; // x coordinate of cable 3, right and down cable
// double y_cable33; // y coordinate of cable 3
// double F_c3; // heat source in cable W/m3


//boundaries temperatures and gradients
// Dirichlet boundary conditions
double Ttop    = 25.0; // temperature at the top of the domain
double Tbottom = 15.0; // temperature at the bottom of the domain
double Tleft   = 15.0; // temperature at the left of the domain
double Tright  = 15.0; // temperature at the right of the domain
double Tini    = 15.0; // initial temperature of the domain
double Tcable  = 90.0; // max temperature of the cable for 

// Neumann boundary conditions
double qtop    = 0.0; // heat flux at the top of the domain
double qbottom = 0.0; // heat flux at the bottom of the domain
double qleft   = 0.0; // heat flux at the left of the domain
double qright  = 0.0; // heat flux at the right of the domain

// kind of boundaries
int b_top    = 1; // 1: Dirichlet, 2: Neumann
int b_bottom = 1; // 1: Dirichlet, 2: Neumann
int b_left   = 1; // 1: Dirichlet, 2: Neumann
int b_right  = 1; // 1: Dirichlet, 2: Neumann

/////////////////////////////////////////
// Parameters to refined grid         //
/////////////////////////////////////////
//  Levels of refined
int n_levels = 6; // number of levels of refinement

//Setting of refinement
// Level 0
double Ax_0 = 4.5 ; // Left top point in x
double Ay_0 = 0.6; // Left top point in y
double Bx_0 = 5.5;   // Right top point in x
double By_0 = 0.6; // Right top point in y
double Cx_0 = 5.5;   // Right bottom point
double Cy_0 = 2.0; // Right bottom point in y
double Dx_0 = 4.5;   // Left bottom point in x
double Dy_0 = 2.0; // Left bottom point in y

// Level 1
double Ax_1 = 4.0 ; // Left top point in x
double Ay_1 = 0.5; // Left top point in y
double Bx_1 = 6.0;   // Right top point in x
double By_1 = 0.5; // Right top point in y
double Cx_1 = 6.0;   // Right bottom point in x
double Cy_1 = 3.0; // Right bottom point in y
double Dx_1 = 4.0;   // Left bottom point in x
double Dy_1 = 3.0; // Left bottom point in y

// Level 2
double Ax_2 = 3.5 ; // Left top point in x
double Ay_2 = 0.4; // Left top point in y
double Bx_2 = 6.5;   // Right top point in x
double By_2 = 0.4; // Right top point in y
double Cx_2 = 6.5;   // Right bottom point in x
double Cy_2 = 4.0; // Right bottom point in y
double Dx_2 = 3.5;   // Left bottom point in x
double Dy_2 = 4.0; // Left bottom point in y

// Level 3
double Ax_3 = 3.0 ; // Left top point in x
double Ay_3 = 0.3; // Left top point in y
double Bx_3 = 7.0;   // Right top point in x
double By_3 = 0.3; // Right top point in y
double Cx_3 = 7.0;   // Right bottom point in x
double Cy_3 = 5.0; // Right bottom point in y
double Dx_3 = 3.0;   // Left bottom point in x
double Dy_3 = 5.0; // Left bottom point in y

// Level 4
double Ax_4 = 2.5 ; // Left top point in x
double Ay_4 = 0.2; // Left top point in y
double Bx_4 = 7.5;   // Right top point in x
double By_4 = 0.2; // Right top point in y
double Cx_4 = 7.5;   // Right bottom point in x
double Cy_4 = 6.0; // Right bottom point in y
double Dx_4 = 2.5;   // Left bottom point in x
double Dy_4 = 6.0; // Left bottom point in y

// Level 5
double Ax_5 = 0.0 ; // Left top point in x
double Ay_5 = 0.0; // Left top point in y
double Bx_5 = DEEP_OF_DOMAIN;   // Right top point in x
double By_5 = 0.0; // Right top point in y
double Cx_5 = DEEP_OF_DOMAIN;   // Right bottom point in x
double Cy_5 = WIDTH_OF_DOMAIN; // Right bottom point in y
double Dx_5 = 0.0;   // Left bottom point in x
double Dy_5 = WIDTH_OF_DOMAIN; // Left bottom point in y

//function to detect the level of refinement of a point give (i,j) coordinates 
// return the steps for the level of refinement
int level(int i, int j) {
    // Define level of refinement based on (x, y) coordinates
    double x =  (double)i * hmin;
    double y =  (double)j * hmin;
    int k = 5; // level of refinement
    if (x >= Ax_0 && x <= Bx_0 && y >= Ay_0 && y <= Dy_0) {k = 0;}
    else if (x >= Ax_1 && x <= Bx_1 && y >= Ay_1 && y <= Dy_1) {k = 1;}
    else if (x >= Ax_2 && x <= Bx_2 && y >= Ay_2 && y <= Dy_2) {k = 2;}
    else if (x >= Ax_3 && x <= Bx_3 && y >= Ay_3 && y <= Dy_3) {k = 3;}
    else if (x >= Ax_4 && x <= Bx_4 && y >= Ay_4 && y <= Dy_4) {k = 4;}
    else if (x >= Ax_5 && x <= Bx_5 && y >= Ay_5 && y <= Dy_5) {k = 5;}
    return pow(2, k);
}

///////////////////////////////////////////////////////////////////////////////
// Function to check if coordinates (x,y) are inside a circle of radio r     //
// and center (x0,y0)                                                       //
///////////////////////////////////////////////////////////////////////////////
int circle(double x, double y, double x0, double y0, double r) {
    return (((x - x0) * (x - x0) + (y - y0) * (y - y0)) <= (r * r));
    }

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Generic function to model the spatial distribution of a "SPATIAL_PARAM"      //
// physic parameter(e.g. density, specific heat, thermal conductivity, etc.) //
///////////////////////////////////////////////////////////////////////////////

double SPATIAL_PARAM(int i, int j, double SP_soil1, double SP_soil2, double SP_xlpe, double SP_cond) {
    // Define SPATIAL PARAMETER based on (x, y) coordinates
        
    double x =  (double)i * hmin;
    double y =  (double)j * hmin;

    double SP_ = SP_soil1; // default value is soil 1 value

    // (x,y) outside soil 2 area - Backfill Soil
    if ((y > (deep_cable + 2*cable_diam)  ||   
        y < (deep_cable - soil2_layer))  || 
        (x > (L/2 +0.4) || x < (L/2 -0.4)))
        {return SP_ ;} // return SP_soil1 if (x,y) is outside soil 2 area
   else {
       SP_ = SP_soil2; //SP_soil is the default value for soil 2

    // Circuit 1, cable 11 , cable 12 and cable 13   - Central circuit
    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        SP_ = SP_cond;
        return SP_; // return SP_cable if (x,y) is inside cable 11, 12 or 13
    }
    if (circle(x, y, x_cable11, y_cable11, cable_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cable_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cable_diam/2)) {
        SP_ = SP_xlpe;
        return SP_; // return SP_xlpe if (x,y) is inside cable 11, 12 or 13
    }
    
    // Circuit 2, cable 21 , cable 22 and cable 23 - Left circuit
    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        SP_ = SP_cond;
        return SP_; // return SP_cable if (x,y) is inside cable 21, 22 or 23
    }
    if (circle(x, y, x_cable21, y_cable21, cable_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cable_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cable_diam/2)) {
        SP_ = SP_xlpe;
        return SP_; // return SP_xlpe if (x,y) is inside cable 21, 22 or 23
    }
    
    // Circuit 3, cable 31 , cable 32 and cable 33 - Right circuit
    // if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
    //     circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
    //     circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
    //     SP_ = SP_cond;
    //     return SP_; // return SP_cable if (x,y) is inside cable 31, 32 or 33
    // }
    // if (circle(x, y, x_cable31, y_cable31, cable_diam/2) || 
    //     circle(x, y, x_cable32, y_cable32, cable_diam/2) || 
    //     circle(x, y, x_cable33, y_cable33, cable_diam/2)) {
    //     SP_ = SP_xlpe;
    //     return SP_; // return SP_xlpe if (x,y) is inside cable 31, 32 or 33
    // }
 
   } // end of else
    return SP_; // return SP_ value
}

///////////////////////////////////////////////////////////////////////////////
// Function to calculate density at a given coordinate (i, j)               //
///////////////////////////////////////////////////////////////////////////////
double rho(int i, int j) {
    // Define density calculation based on (x, y) coordinates
    return SPATIAL_PARAM(i, j, rho_soil1, rho_soil2, rho_xlpe, rho_cond);
}

///////////////////////////////////////////////////////////////////////////////
// Function to calculate specific heat Ce at a given coordinate (i, j)      //
///////////////////////////////////////////////////////////////////////////////
double Ce(int i, int j) {
    // Define Specific Heat calculation based on (x, y) coordinates
    return SPATIAL_PARAM(i, j, Ce_soil1, Ce_soil2, Ce_xlpe, Ce_cond);;
}
///////////////////////////////////////////////////////////////////////////////
// Function to calculate thermal conductivity "K" at a coordinate (i, j)     //
///////////////////////////////////////////////////////////////////////////////
double K(int i, int j) {
    // Define thermal conductivity calculation based on (x, y) coordinates
    return SPATIAL_PARAM(i, j, k_soil1, k_soil2, k_xlpe, k_cond);;
}
///////////////////////////////////////////////////////////////////////////////
// Function to calculate heat source "F" at a given coordinate (i, j)          //
///////////////////////////////////////////////////////////////////////////////
double F(int i, int j, int step, double d_t) {
    // Define source heat calculation based on (x, y) coordinates
    double F_ = 0.0; // Assumed that only conductors are generating heat
    
    double x =  (double)i * hmin;
    double y =  (double)j * hmin;

    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        return F_c1*d_t/(Ce(i,j) * rho(i,j)); // return F_cable if (x,y) is inside cable 11, 12 or 13
    }

    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
         circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
         circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
         return F_c2*d_t/(Ce(i,j) * rho(i,j)); // return F_cable if (x,y) is inside cable 21, 22 or 23
     }
    
    // if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
    //     circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
    //     circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
    //     return F_c3 *d_t/(Ce(i,j) * rho(i,j))
    //}
    return F_;
}

///////////////////////////////////////////////////////////////////////////////
// Function to return initial temperature given coordinate (x, y)           //
///////////////////////////////////////////////////////////////////////////////
double U(int i, int j) {
    double x =  (double)i * hmin;
    double y =  (double)j * hmin;

    // Define heat source calculation based on (x, y) coordinates
    return Tini; // Initial temperature
}  

///////////////////////////////////////////////////////////////////////////////
// Function to calculate the probabilities of the walker to move to the      //
// next position                                                            //
/////////////////////////////////////////////////////////////////////////////// 
int calculate_probabilities( int i, int j, int steps, double d_t) {
    // Define the probabilities of the walker to move to the next position
    // based on the thermal conductivity of the material at the next position
    double Ce_E, Ce_W, Ce_S, Ce_N , Ce_P;    // Coefficients of the temperature terms of the heat equation at interior node
    double p_P, p_E, p_W, p_N, p_S;      // Probabilities of the walker to move to the next position
    double K_E, K_W, K_N, K_S, K_P; // thermal conductivity around the point P --> (i,j)
    double F_e, F_w, F_n, F_s, F_p; // linearly interpolated values of thermal conductivity
    double F_c ; // linearly interpolated values of thermal conductivity of central point
    double rho_E, rho_W, rho_N, rho_S, rho_P; // linearly interpolated values of thermal conductivity
    double omega = 1.0*d_t / (2.0*hmin*steps); //1.0*d_t / (hmin*hmin*steps*steps); // time step divided by h^2

    int verbose = VERBOSE3; // 1 = print probabilities, 0 = don't print probabilities

    // Values of thermal conductivity 
    K_P = K(i  ,j);    // thermal conductivity at the point P --> (i,j)
    K_W = K(i-steps,j);    // thermal conductivity at the point W --> (i-1,j)  
    K_E = K(i+steps,j);    // thermal conductivity at the point E --> (i+1,j)
    K_S = K(i  ,j+steps ); // thermal conductivity at the point S --> (i,j+1)
    K_N = K(i  ,j-steps);  // thermal conductivity at the point N --> (i,j-1)

    // Parametros en el punto P
    rho_P = rho(i  , j) ;    // densidad en el punto P --> (i,j)
    Ce_P = Ce(i  ,j);    // calor especifico en el punto P --> (i,j)
    // verificar silos valores de K son iguales 
    if (K_P == K_E && K_P == K_W && K_P == K_N && K_P == K_S) {
        F_c = K_P / (rho_P * Ce_P) * omega;
        p_P = (1.0 - 4.0*F_c);
        p_E = (1.0 - (3.0*F_c)); // probability to move to the left
        p_W = (1.0 - (2.0*F_c)); // probability to move to the right
        p_N = (1.0 - F_c); // probability to move to the top
        p_S = 1.0; // probability to move to the bottom
    } // else
    else {
        // Valores de rho 
        rho_W = rho(i-steps, j); // densidad en el punto W --> (i-1,j)
        rho_E = rho(i+steps, j); // densidad en el punto E --> (i+1,j)
        rho_S = rho(i  , j+steps); // densidad en el punto S --> (i,j+1)
        rho_N = rho(i  , j-steps); // densidad en el punto N --> (i,j-1)

        // Valores de Ce 
        Ce_W = Ce(i-steps, j); // calor especifico en el punto W --> (i-1,j)
        Ce_E = Ce(i+steps, j); // calor especifico en el punto E --> (i+1,j)
        Ce_S = Ce(i  , j+steps); // calor especifico en el punto S --> (i,j+1)
        Ce_N = Ce(i  , j-steps); // calor especifico en el punto N --> (i,j-1)

        // Valores de F (coeficientes) interpolados linealmente
        F_w = 2* (K_W + K_P) /((rho_W + rho_P)*(Ce_W + Ce_P)) * omega ; // valor de F entre W y P
        F_e = 2* (K_E + K_P) /((rho_E + rho_P)*(Ce_E + Ce_P)) * omega ; // valor de F entre E y P
        F_n = 2* (K_N + K_P) /((rho_N + rho_P)*(Ce_N + Ce_P)) * omega ; // valor de F entre N y P
        F_s = 2* (K_S + K_P) /((rho_S + rho_P)*(Ce_S + Ce_P)) * omega ; // valor de F entre S y P
        F_p = (K_E + K_W + K_N + K_S) /((rho_P* Ce_P)) * omega ;        // valor de F en P

        // Cumulative probabilities to next position transition
        p_P = (1 - F_p);
        p_E = (1 - (F_p - F_e)); // probability to move to the left
        p_W = (1 - (F_p - F_e - F_w)); // probability to move to the right
        p_N = (1 - (F_p - F_e - F_w - F_n)) ;// probability to move to the top
        p_S = 1.00; // probability to move to the bottom
    }
    
    // print p_E, p_W, p_N and p_S conditioned to verbose in a single row
    if (verbose == 1) {
        fprintf(stderr," %d: Probabilities:  %.5f %.5f %.5f %.5f %.5f\n ", omp_get_thread_num(), p_P, p_E, p_W, p_N, p_S);
    }

    // random number between 0 and 1 uniformly distributed using rand() function (RANDOM_LIB = 1) or randlib library (RANDOM_LIB = 0)
    double r; // random number
    if (RANDOM_LIB == 1) { 
        r = (double)rand() / (double)RAND_MAX;}
    else { 
        // from the interval [0, 1) using the ranlib library
        r =  genunf ( 0, 1); }
    
    // double r =  genunf ( 0, 1); 
    if (verbose == 1) {fprintf(stderr," random number: %.6f\n", r);}
    // random directions or walker
    if (r <= p_P) { // No move 
        return 0;
    }
    if (r <= p_E) { // move to the right
        return 1;
    }
    if (r <= p_W) { // move to the left
        return 2;
    }
    if (r <= p_N) { // move to the top
        return 3;
    }
    return 4; // move to the bottom
}

//////////////////////////////////////////////////////////////////////////////////////////
// Functions to calculate the temperature of a single walker in the domain              //
//////////////////////////////////////////////////////////////////////////////////////////
double single_walk(int start_i, int start_j) {
    // single random walk
    int i = start_i; // initial walker position
    int j = start_j; // initial walker position
    int direction; // random direction
    double temp_walker = 0 ; // scored walker temperature
    int no_border = 1; // 0 = border touched, 1 = no border touched
    int verbose = VERBOSE2 ; // 1 = print walker position, 0 = don't print walker position
    int walk_steps=0; // counter for the walk steps number 
    int steps; // number of steps in each walk accordance with the refinement level
    double time = 0.0; // time point for the walker
    double d_t = d_t_0; // time step
    double tmax = TMAX; // maximum time for the walker
    
    // loop over the walk length 
    while (no_border == 1 && time <= tmax) {
        steps = level(i,j); // steps in the level of refinement
        d_t = (steps) * d_t_0; // time step size for cooper conductivity
        time += d_t; // update the time
        
        walk_steps++; // update the walk steps number
        //Add heat source to the temperature
        temp_walker += F(i,j, steps, d_t)  ; // *steps;
        //Calculate the probabilities and move the walker 
        direction = calculate_probabilities(i, j, steps, d_t);
        switch (direction) {
            case 0:
                //no_border = 0; // No border
                //temp_walker += U(i,j);
                break; // No move
            case 1:
                i = i + steps; // move to the right
                if (i > nx - steps) { // Right border
                    if (verbose == 1) {fprintf(stderr," %d: Right Border touched! \n ", omp_get_thread_num());}
                    if (b_right == 1) { // Dirichlet border
                        temp_walker += Tright; 
                        no_border = 0; // No border
                    } else { // Neumann border (with q=0)
                        i= nx-steps; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 2:
                i = i - steps; // move to the left
                if (i < steps) {  // Left border
                    if (verbose == 1) {fprintf(stderr," %d: Left Border touched! \n ", omp_get_thread_num());}
                    if (b_left == 1) { // Dirichlet border
                        temp_walker += Tleft;
                        no_border = 0; // No border 
                    } else { // Neumann border  (with q=0)
                        i =steps; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    }
                }
                break;
            case 3:
                j = j - steps; // move to the top
                if (j < steps) { // Top border
                    if (verbose == 1) {fprintf(stderr," %d: Top Border touched! \n ", omp_get_thread_num());}
                    if (b_top == 1) { // Dirichlet border
                        temp_walker += Ttop;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= steps; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 4:
                j = j + steps; // move to the bottom
                if (j > ny - steps) { // Bottom border
                    if (verbose == 1) {fprintf(stderr," %d: Bottom Border touched! \n ", omp_get_thread_num());}
                    if (b_bottom == 1) { // Dirichlet border
                        temp_walker += Tbottom;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= ny-steps; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break; // end of the switch
        }

    // print if no_border == 0 conditioned to verbose
/*
    if (VERBOSE4) {
        fprintf(stderr,"num steps: %d \n", steps);
        // print time, i, j, temp_walker, walk_steps
        fprintf(stderr,"%d: Time: %.2f, i: %d, j: %d, Temp: %.2f, Walk steps: %d\n", omp_get_thread_num(), time, i, j, temp_walker, walk_steps);
        if (no_border == 0 ) {
            fprintf(stderr,"%d: Walker out of bounds!\n", omp_get_thread_num());
        }
    }
*/
    } // end of the walk loop

    if (no_border == 1) { // if the walker is still inside the domain at the end of the walk
        temp_walker += U(i,j); 
       }
    if (VERBOSE4) {
        fprintf(stderr,"%d: Time: %.2f, i: %d, j: %d, Temp: %.2f, Walk steps: %d\n", omp_get_thread_num(), time, i, j, temp_walker, walk_steps);
        fprintf(stderr,"%d: Walker scored Temp: %.2f\n", omp_get_thread_num(), temp_walker);
        fprintf(stderr,"%d: Walker walk steps: %d\n", omp_get_thread_num(), walk_steps);
        }
    return temp_walker; // return the calculated temperature of the walker
}

//////////////////////////////////////////////////////////////////////////////////////////
// Main function                                                                       //
//////////////////////////////////////////////////////////////////////////////////////////
int main(void)
{
    double start = omp_get_wtime(); // start counting of time
    //////////////////////////////////////////////
    // Initialize the random numbers generators //
    //////////////////////////////////////////////
    // Initialize the random number generator rand() with the current time
    srand(time(NULL));

    unsigned long HAUSNUMERO = 12345;
    char phrase[] = "randomizer";
    char  buf[] = "randomizer1";
    int seed1;
    int seed2;
 
    initialize ( ); // initialize the random number generators
    // use sprintf() function to add string "phrase" to time and convert to a string in buf
    sprintf ( buf, "%d%s", time(NULL), phrase );
    // Set the seeds based on the phrase.
    phrtsd ( buf, &seed1, &seed2 );
    // Initialize all generators.
    set_initial_seed ( seed1, seed2 );
    
    ///////////////////////////////////////////////
    //           space step size h               //
    //////////////////////////////////////////////
    hmin= L / ( (double) nx - 1); // nx-1 because we have nx points in the domain and nx-1 intervals

    ///////////////////////////////////////////////
    //           time step size d_t              //
    //////////////////////////////////////////////
    d_t_0 = Factor_dt * (2* hmin ) * rho_cond * Ce_cond / (4 * k_cond); // time step size for cooper conductors

    //////////////////////////////////////////////
    // Cable geometry and arrangement  parameters//
    //////////////////////////////////////////////
    double heigh_central_cable; // heigh of the central cable
    double separat_lateral_cables; // separation between lateral cables
    int type_cbl_arrange= TYPE_CABLE_ARRANGEMENT; // type of cable arrangement

    // Horizontal Cables Arrangement
    if (type_cbl_arrange == 1){ // Horizontal Cables Arrangement
        heigh_central_cable = 0.0;
        separat_lateral_cables= cable_diam +0.5;
    }
    // Trifoil Cable Arrangement
    if (type_cbl_arrange == 2){ // Trifoil Cable Arrangement
        heigh_central_cable = cable_diam*(0.866025) + 0.001;
        separat_lateral_cables= 0.001;
    }
    // print type of cable arrangement
    if (VERBOSE) {
        fprintf(stderr,"Cable arrangement type: %d\n", type_cbl_arrange);
    }

    // Initialize values for circuits
    //Circuit 1 - three phase arrangement in trifoil - Central circuit
    x_cable11 = L/2 + separation_circuits/2 ; // x coordinate of cable 1 , central and upper cable
    y_cable11 = deep_cable - heigh_central_cable - cable_diam; // y coordinate of cable 1
    x_cable12 = x_cable11 - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    y_cable12 = deep_cable - cable_diam; // y coordinate of cable 2
    x_cable13 = x_cable11 + cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    y_cable13 = deep_cable - cable_diam; // y coordinate of cable 3
        
    // current in circuit 1
    double Ic1 = 328.6; // current in cable
    int icable1= (int) (x_cable11/hmin); // cable number
    int jcable1= (int) (y_cable11/hmin); // cable number
    int i = (int) (x_cable11 / hmin);
    int j = (int) (y_cable11 / hmin);
    
    F_c1 = 0.0 ;// (1.0 * Ic1 * Ic1 * res_cond/(3.1416*cond_diam*cond_diam/4)) ; // heat source in cable W/m3 --> 1.05 is for insulation losses
    //F_c1=(1.05 * Ic1 * Ic1 * res_cond/(3.1416*cond_diam*cond_diam/4)) ; // heat source in cable W/m3 --> 1.05 is for insulation losses
    //Circuit 2 - three phase arrangement in trifoil - Left circuit
    x_cable21 = x_cable11 - separation_circuits; // x coordinate of cable 1 , central and upper cable
    y_cable21 = y_cable11; // y coordinate of cable 1
    x_cable22 = x_cable21 - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    y_cable22 = deep_cable - cable_diam; // y coordinate of cable 2
    x_cable23 = x_cable21 +  cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    y_cable23 = deep_cable - cable_diam; // y coordinate of cable 3
    
    // current in circuit 2
    double Ic2 = 490.0; //690.0; // 328.6; // current in cable
    int icable2= (int) (x_cable21/hmin); // cable number
    int jcable2= (int) (y_cable21/hmin); // cable number
    F_c2 = 1.05* Ic2 * Ic2 * res_cond/(3.1416*cond_diam*cond_diam/4); // heat source in cable W/m3 --> 1.05 is for insulation losses
    //F_c2 = 1.05* Ic2 * Ic2 * res_cond/(3.1416*cond_diam*cond_diam/4) ; // heat source in cable W/m3 --> 1.05 is for insulation losses 

    //Circuit 3 - three phase arrangement in trifoil - Right circuit
    // x_cable31 = x_cable11 + separation_circuits; // x coordinate of cable 1 , central and upper cable
    // y_cable31 = y_cable11; // y coordinate of cable 1
    // x_cable32 = x_cable31  - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    // y_cable32 = deep_cable - cable_diam/2; // y coordinate of cable 2
    // x_cable33 = x_cable31  + cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    // y_cable33 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // // current in circuit 3
    // double Ic3 = 0.0; // current in cable
    // int icable3= (int) (x_cable31/hmin); // cable number
    // int jcable3= (int) (y_cable31/hmin); // cable number
    // F_c3 = 1.05*Ic3*Ic3*res_cond/(cable_area); // heat source in cable W/m3

    // Point to calculate the temperature
    ///////////////////////////////////////
    ///////////////////////////////////////
    double x_point = x_cable21;          //   L/2; // x coordinate of the point
    double y_point = y_cable21 ; // - cable_diam/2-0.002; // y coordinate of the point
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Initial point of walkers
    int start_i = (int) (x_point/hmin); // starting x coordinate
    int start_j = (int) (y_point/hmin); // starting y coordinate
    int steps = level(start_i, start_j); // steps in the level of refinement

    // loops
    int verbose = VERBOSE ; // 1 = print walker position, 0 = don't print walker position
    int n_loops=0;
    double T_walkers=0.0;
    double sum_temperature = 0.0;
    double Sum_Temp = 0.0;
    double sum_squared_temperature = 0.0;
    double Sum_Sqr_Temp = 0.0;
    double average_temperature=0.0;
    double variance = 0.0;
    
    // set the number of threads
    int n_walkers = (N_WALKERS == 0) ? omp_get_max_threads() : N_WALKERS; // number of walkers
    omp_set_num_threads(n_walkers);

    // print a descriptive wording about the simulation including some important data about it in the screen
    printf("Calculation of thermal heat of electric conductors in a cable arrangement\n");
    if (type_cbl_arrange == 1) {
        printf("Cable arrangement: 3 cables in a row\n");
    } else if (type_cbl_arrange == 2) {
        printf("Cable arrangement: 3 cables in a trifoil\n");
    }
    printf("Number of parallel walkers: %d\n", n_walkers);
    printf("Number of loops: %d\n", MAX_NUMBER_OF_LOOPS);
    // print starting point of walkers
    printf("Starting point of walkers: (%d,%d)\n", start_i, start_j);
    // print h and L
    printf("hmin: %.6f\n", hmin);
    printf("hmax: %.6f\n", pow(2,(n_levels-1))*hmin);
    printf("L: %.6f\n", L);

    ////////////////////////////////////////////////
    // loop over the walk length for multi walkers//
    ////////////////////////////////////////////////
    
    //advance control
    double VarToAvrg= 99999.9;
    double percentage_advance = 0; // percentage of advance
    double perc_advance_step = 1.0; // percentage step to print advance
    int advance_step = (int) MAX_NUMBER_OF_LOOPS* (perc_advance_step/100.0); // percentage to print advance
    int print_advance = ((int)(advance_step * n_walkers * nx) > 1000) ? PRINT_ADVANCE : 0 ; //dont print for small calcs
    if (print_advance) {
    printf("percentage of advance: %0.2f%%  - T_avg: %.6f - variance-to-mean ratio: %.6f \r", 
                    percentage_advance, average_temperature, VarToAvrg);}
    
    // loop over the walk length for multi walkers
    while ((n_loops < MAX_NUMBER_OF_LOOPS)
           && ((VarToAvrg > MAX_VARIANCE_TO_MEAN_RATIO)
           || (n_loops < 4))) {

        // while (n_loops < MAX_NUMBER_OF_LOOPS) {
        n_loops++;

        // start #pragma parallel region
        sum_temperature = 0.0;
        sum_squared_temperature = 0.0;

         #pragma omp parallel reduction(+:sum_temperature) \
                 reduction(+:sum_squared_temperature)      \
                 private (T_walkers) 
        {      
            #pragma omp for nowait //schedule(dynamic)
                for(int k = 0; k < n_walkers; k++) {
                    T_walkers = single_walk(start_i, start_j);
                    sum_temperature += T_walkers;
                    sum_squared_temperature += T_walkers*T_walkers;
                    } // close #pragma omp for
        }
        #pragma omp barrier
        // Calculate variance and average temperature
        Sum_Temp += sum_temperature;
        Sum_Sqr_Temp += sum_squared_temperature;

        average_temperature = (Sum_Temp) / ((double) n_walkers * (double) n_loops);
        variance = (Sum_Sqr_Temp / (n_walkers*n_loops) - (average_temperature * average_temperature))/(n_walkers * n_loops);
        VarToAvrg = (variance/average_temperature);
        if (print_advance) { // print advance percentage every 5% of MAX_NUMBER_OF_LOOPS 
            if (n_loops % advance_step == 0) {
                percentage_advance = (double) ((n_loops*100.0)/MAX_NUMBER_OF_LOOPS);
                average_temperature = (Sum_Temp) / (n_walkers * n_loops);
                printf("percentage of advance: %.2f%%  - T_avg: %.6f - variance-to-mean ratio: %.6f \r", 
                    percentage_advance, average_temperature, VarToAvrg); 
            }
        } // close if
    
    } // close while loop

    // Calculate average temperature
 
    double end = omp_get_wtime();
    // print time of execution
    printf("\n time: %f\n",end-start);
    printf("Average Temperature: %.8f\n", average_temperature);
    printf((" at time: %.8f\n"), TMAX);
    printf((" delta t_0: %.8f\n"), d_t_0);
    printf("Standard Deviation: %.8f\n", sqrt(variance));
    printf("variance-to-mean ratio: %.6f\n", VarToAvrg);
    printf("Walkers: %d \n", (n_walkers*n_loops));
    printf("Loops: %d \n", (n_loops));
    printf("parallel walkers: %d \n", (n_walkers));

    // write in a single text row Average Temperature, standard deviation and  n_walkers*n_loops to file to a file (new or add to existent) 
    FILE *fp; // file pointer
    fp = fopen("results_ex3_transient.txt", "a"); // open file to append
    fprintf(fp, "%.8f   %.8f   %d    %.8f %.8f  %d \n", average_temperature, sqrt(variance), n_walkers*n_loops, hmin, d_t_0 ,nx) ;
    fclose(fp);
   return 0;
}
