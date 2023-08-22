// Heat equation in 2d with random numbers for steady state
//  gcc -o heat_random_v1.exe heat_random_v1.c  ranlib.c rnglib.c -fopenmp
// heat_random_v1.exe
// Main parameters
#define N_WALKERS 0 // number of parallel walkers 0 for maximum hardware number of threads
#define MAX_VARIANCE_TO_MEAN_RATIO 0.00075 // maximum variance to mean ratio to stop the simulation
#define DEEP_OF_DOMAIN 2.0 // deep of domain in meters
#define WIDTH_OF_DOMAIN 2.0 // width of domain in meters
#define NUMBER_OF_POINTS_IN_X 10000
#define NUMBER_OF_POINTS_IN_Y 10000
#define MAX_NUMBER_OF_LOOPS 1000// maximum number of loops
#define DEEP_CABLE 1.0 // deep of cable in meters
#define SEP_CABLES 0.2  // separation between cables in meters between three phase circuits
#define TYPE_CABLE_ARRANGEMENT 2 // 1 for horizontal arrangement, 2 for trifoil arrangement
#define BACKFILL_HEIGHT_LAYER 0.25 // height of backfill soil layer in m
#define VERBOSE 0 // 1 to print detailed results of total results in the screen
#define VERBOSE2 0 // 1 to print detailed results of total results in the screen
#define VERBOSE3 0 // 1 to print detailed results of probabilities in the screen
#define VERBOSE4 0 // 1 to print detailed results of heat sources in the screen
#define PRINT_ADVANCE 1 // 1 to print advance of the simulation in the screen
#define RANDOM_LIB 1 // 1 to use native random library, 0 to use ranlib library

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

double h; // spatial step
int n_steps; // number of steps in each walk

// Aluminum 630 mm2 data for calculations
double cable_area= 0.000630;// cross area of cable 630 mm2 in squered meters
double cable_diam = 0.0552; // diameter of cable 630 mm2 - fuente: http://es.lemeicable.com/product/108.html
double cond_diam =  0.0301; // diameter of conductor 630 mm2 - fuente: http://es.lemeicable.com/product/108.html 
double res_cond = 0.000064; // ohm/m  - electric resistivity of cable 630 mm2 0.0469 ohm/km at 20°C DC - fuente: http://es.lemeicable.com/product/108.html
double rho_cond = 2700.0; // density of aluminum cable in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_cond = 900.0; // specific heat of aluminum cable in J/(kg K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double k_cond = 237.0; // thermal conductivity of aluminum cable in W/(m.K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html

// datos XLPE
double rho_xlpe = 920.0; // density of XLPE in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_xlpe = 2174.0; // specific heat of XLPE in J/(kg K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
double k_xlpe = 0.28; // thermal conductivity of XLPE in W/(m.K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
double F_xlpe = 0.0; // heat source in cable W/m3

// Native Soil thermal data for calculations
double rho_soil1 = 2000.0; // density of soil_1 in kg/m3 - fuente: https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil1 = 800.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil1 = 0.125; // thermal conductivity of soil 1 in W/(m.K) (resistivity = 8)- fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5
double F_soil1 = 0.0; // heat source in soil 1 W/m3

// Backfill Soil thermal data for calculations
double rho_soil2 = 1500.0; // density of soil_1 in kg/m3 - fuente:  https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil2 = 900.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil2 = 0.333; // thermal resistivity of soil 1 in W/(m.K) (resistivity = 3) - fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5 
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

//Circuit 3 - three phase arrangement in trifoil - Right circuit
double x_cable31; // x coordinate of cable 1 , central and upper cable
double y_cable31; // y coordinate of cable 1
double x_cable32; // x coordinate of cable 2, left and down cable
double y_cable32; // y coordinate of cable 2
double x_cable33; // x coordinate of cable 3, right and down cable
double y_cable33; // y coordinate of cable 3
double F_c3; // heat source in cable W/m3


//boundaries temperatures and gradients
// Dirichlet boundary conditions
double Ttop    = 20.0; // temperature at the top of the domain
double Tbottom = 15.0; // temperature at the bottom of the domain
double Tleft   = 20.0; // temperature at the left of the domain
double Tright  = 20.0; // temperature at the right of the domain
double Tini    = 20.0; // initial temperature of the domain
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
        
    double x =  (double)i * h;
    double y =  (double)j * h;

    double SP_ = SP_soil1; // default value is soil 1 value

    // (x,y) outside soil 2 area - Backfill Soil
    if ((y > (deep_cable + cable_diam)  ||   
        y < (deep_cable - soil2_layer))  || 
        (x > (3*L/4) || x < (L/4)))
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
    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        SP_ = SP_cond;
        return SP_; // return SP_cable if (x,y) is inside cable 31, 32 or 33
    }
    if (circle(x, y, x_cable31, y_cable31, cable_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cable_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cable_diam/2)) {
        SP_ = SP_xlpe;
        return SP_; // return SP_xlpe if (x,y) is inside cable 31, 32 or 33
    }
 
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
double F(int i, int j) {
    // Define source heat calculation based on (x, y) coordinates
    double F_ = 0.0; // Assumed that only conductors are generating heat
    
    double x =  (double)i * h;
    double y =  (double)j * h;

    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        F_ = F_c1;
        return F_;
    }

    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        F_ =  F_c2;
        return F_;
    }
    
    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        F_ =  F_c3;
        return F_;
    }
    return F_;
}

///////////////////////////////////////////////////////////////////////////////
// Function to return initial temperature given coordinate (x, y)           //
///////////////////////////////////////////////////////////////////////////////
double U(int i, int j) {
    double x =  (double)i * h;
    double y =  (double)j * h;

    // Define heat source calculation based on (x, y) coordinates
    return Tini; // Initial temperature
}  

///////////////////////////////////////////////////////////////////////////////
// Function to calculate the probabilities of the walker to move to the      //
// next position                                                            //
/////////////////////////////////////////////////////////////////////////////// 
int calculate_probabilities( int i, int j) {
    // Define the probabilities of the walker to move to the next position
    // based on the thermal conductivity of the material at the next position
    double p_e, p_ew, p_ewn , p_ewns ;

    int verbose = VERBOSE3; // 1 = print probabilities, 0 = don't print probabilities
        
    double Kxy =   K(i     , j);
    double Kxm_1y = K(i - 1 , j);
    double Kxp_1y = K(i + 1 , j);
    double Kxym_1 = K(i     , j - 1 );
    double Kxyp_1 = K(i     , j + 1);
    
    // if all "Ks" are the same then all probabilities are equal to 0.25
    if (Kxy == Kxm_1y && Kxy == Kxp_1y && Kxy == Kxym_1 && Kxy == Kxyp_1) {
        p_e  = 0.25;
        p_ew = 0.5;
        p_ewn = 0.75;
        p_ewns = 1;
    }
    else {
        double halfKxm1yKxy = 0.5 * (Kxm_1y + Kxy);
        double halfKxp1yKxy = 0.5 * (Kxp_1y + Kxy);

        double halfKxym1Kxy = 0.5 * (Kxym_1 + Kxy);
        double halfKxyp1Kxy = 0.5 * (Kxyp_1 + Kxy);

        p_e  = 0.25 + ( halfKxp1yKxy - halfKxm1yKxy)/(8*Kxy);

        p_ew = 0.5;
        p_ewn = 0.75 +(-halfKxyp1Kxy + halfKxym1Kxy)/(8*Kxy); //corregido -ojo
        p_ewns = 1;
    }
    // print p_e, P_ew, P_ewn, P_ewns and the sum of all probabilities
    if (verbose == 1) { 
        fprintf(stderr,"p_e = %.6f, p_ew = %.6f, p_ewn = %.6f, p_ewns = %.6f\n", p_e, p_ew, p_ewn, p_ewns);
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
    if (r <= p_e) { // move to the right
        return 1;
    }
    if (r <= p_ew) { // move to the left
        return 2;
    }
    if (r <= p_ewn) { // move to the top
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
    double temp_walker = 0.0 ; // scored walker temperature
    double Sf = 0.0 ; // Sum of heat source score
    int no_border = 1; // 0 = border touched, 1 = no border touched
    int verbose = VERBOSE2 ; // 1 = print walker position, 0 = don't print walker position
    int walk_steps=0; // counter for the walk steps number 

    // loop over the walk length 
    while (no_border == 1) {
        walk_steps++;
        //Add heat source to the temperature
        Sf += F(i,j); 
        //Calculate the probabilities and move the walker 
        direction = calculate_probabilities(i, j);
        switch (direction) {
            case 1:
                i++; // move to the right
                if (i > nx - 1) { // Right border
                    if (verbose == 1) {fprintf(stderr," %d: Right Border touched! \n ", omp_get_thread_num());}
                    if (b_right == 1) { // Dirichlet border
                        temp_walker += Tright + Sf; 
                        no_border = 0; // No border
                    } else { // Neumann border (with q=0)
                        i= nx-1; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 2:
                i--; // move to the left
                if (i < 2) {  // Left border
                    if (verbose == 1) {fprintf(stderr," %d: Left Border touched! \n ", omp_get_thread_num());}
                    if (b_left == 1) { // Dirichlet border
                        temp_walker += Tleft + Sf;
                        no_border = 0; // No border 
                    } else { // Neumann border  (with q=0)
                        i =2; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    }
                }
                break;
            case 3:
                j--; // move to the top
                if (j < 2) { // Top border
                    if (verbose == 1) {fprintf(stderr," %d: Top Border touched! \n ", omp_get_thread_num());}
                    if (b_top == 1) { // Dirichlet border
                        temp_walker += Ttop + Sf;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= 2; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 4:
                j++; // move to the bottom
                if (j > ny - 1) { // Bottom border
                    if (verbose == 1) {fprintf(stderr," %d: Bottom Border touched! \n ", omp_get_thread_num());}
                    if (b_bottom == 1) { // Dirichlet border
                        temp_walker += Tbottom + Sf;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= ny-1; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
        }
    // print if no_border == 0 conditioned to verbose
    if (no_border == 0 && VERBOSE4) {
        fprintf(stderr,"%d: Walker out of bounds!\n", omp_get_thread_num());}
    } // end of the walk loop

    if (no_border == 1) { // if the walker is still inside the domain at the end of the walk
        temp_walker += U(i,j) + Sf; 
       }
    if (VERBOSE4) {
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
    h= (double) L/( (double)nx - 1); // nx-1 because we have nx points in the domain and nx-1 intervals

    //////////////////////////////////////////////
    // Cable geometry and arrangement  parameters//
    //////////////////////////////////////////////
    double heigh_central_cable; // heigh of the central cable
    double separat_lateral_cables; // separation between lateral cables
    int type_cbl_arrange= TYPE_CABLE_ARRANGEMENT; // type of cable arrangement

    // Horizontal Cables Arrangement
    if (type_cbl_arrange == 1){ // Horizontal Cables Arrangement
        heigh_central_cable = 0.0;
        separat_lateral_cables= cable_diam;
    }
    // Trifoil Cable Arrangement
    if (type_cbl_arrange == 2){ // Trifoil Cable Arrangement
        heigh_central_cable = cable_diam*(0.866025);
        separat_lateral_cables= 0.0;
    }
    // print type of cable arrangement
    if (VERBOSE) {
        fprintf(stderr,"Cable arrangement type: %d\n", type_cbl_arrange);
    }

    // Initialize values for circuits
    //Circuit 1 - three phase arrangement in trifoil - Central circuit
    x_cable11 = L/2; // x coordinate of cable 1 , central and upper cable
    y_cable11 = deep_cable - heigh_central_cable -0.5*cable_diam; // y coordinate of cable 1
    x_cable12 = x_cable11 - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    y_cable12 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable13 = x_cable11 + cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    y_cable13 = deep_cable - cable_diam/2; // y coordinate of cable 3
        
    // current in circuit 1
    double Ic1 = 330; // current in cable
    int icable1= (int) (x_cable11/h); // cable number
    int jcable1= (int) (y_cable11/h); // cable number
    F_c1 = 1.15* Ic1*Ic1*res_cond/(cable_area) *h*h /(4.0*K(icable1,jcable1)); // heat source in cable W/m3 --> 1.15 is for insulation losses
   
    //Circuit 2 - three phase arrangement in trifoil - Left circuit
    x_cable21 = x_cable11 - separation_circuits; // x coordinate of cable 1 , central and upper cable
    y_cable21 = y_cable11; // y coordinate of cable 1
    x_cable22 = x_cable21 - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    y_cable22 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable23 = x_cable21 +  cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    y_cable23 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // current in circuit 2
    double Ic2 = 330; // current in cable
    int icable2= (int) (x_cable21/h); // cable number
    int jcable2= (int) (y_cable21/h); // cable number
    F_c2 = 1.15*Ic2*Ic2*res_cond/(cable_area) *h*h /(4.0*K(icable2,jcable2)); // heat source in cable W/m3 --> 1.15 is for insulation losses

    //Circuit 3 - three phase arrangement in trifoil - Right circuit
    x_cable31 = x_cable11 + separation_circuits; // x coordinate of cable 1 , central and upper cable
    y_cable31 = y_cable11; // y coordinate of cable 1
    x_cable32 = x_cable31  - cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 2, left and down cable
    y_cable32 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable33 = x_cable31  + cable_diam*(0.5 + separat_lateral_cables); // x coordinate of cable 3, right and down cable
    y_cable33 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // current in circuit 3
    double Ic3 = 330.0; // current in cable
    int icable3= (int) (x_cable31/h); // cable number
    int jcable3= (int) (y_cable31/h); // cable number
    F_c3 = 1.15*Ic3*Ic3*res_cond/(cable_area) *h*h /(4.0*K(icable3,jcable3)); // heat source in cable W/m3 --> 1.15 is for insulation losses 

    // Point to calculate the temperature
    ///////////////////////////////////////
    ///////////////////////////////////////
    double x_point = x_cable11;          //   L/2; // x coordinate of the point
    double y_point = y_cable11 - 0.5*(cable_diam - cond_diam );          //  deep_cable - 1.15*cable_diam; // y coordinate of the point
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Initial point of walkers
    int start_i = (int) (x_point/h); // starting x coordinate
    int start_j = (int) (y_point/h); // starting y coordinate
    
    // loops
    int verbose = VERBOSE ; // 1 = print walker position, 0 = don't print walker position
    int n_loops=0;
    int T_walkers=0.0;
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

    ////////////////////////////////////////////////
    // loop over the walk length for multi walkers//
    ////////////////////////////////////////////////
    
    //advance control
    double VarToAvrg= 99999.9;
    int percentage_advance = 0; // percentage of advance
    int perc_advance_step = 1; // percentage step to print advance
    int advance_step = MAX_NUMBER_OF_LOOPS* perc_advance_step/100; // percentage to print advance
    int print_advance = ((advance_step * n_walkers * nx) > 1000) ? PRINT_ADVANCE : 0 ; //dont print for small calcs
    if (print_advance) {
    printf("percentage of advance: %d%%  - T_avg: %.6f - variance-to-mean ratio: %.6f \r", 
                    percentage_advance, average_temperature, VarToAvrg);}
    
    // loop over the walk length for multi walkers
    while ((n_loops < MAX_NUMBER_OF_LOOPS)
           && ((VarToAvrg > MAX_VARIANCE_TO_MEAN_RATIO)
           || (n_loops < advance_step))) {
        n_loops++;

        // start #pragma parallel region
        sum_temperature = 0.0;
        sum_squared_temperature = 0.0;

        #pragma omp parallel for reduction(+:sum_temperature) \
                 reduction(+:sum_squared_temperature)      
        for(int k = 0; k < n_walkers; k++) {
            T_walkers = single_walk(start_i, start_j);
            sum_temperature += T_walkers;
            sum_squared_temperature += T_walkers*T_walkers;
            } // close #pragma omp for
        
        // Calculate variance and average temperature
        Sum_Temp += sum_temperature;
        Sum_Sqr_Temp += sum_squared_temperature;

        average_temperature = (Sum_Temp) / (n_walkers * n_loops);
        variance = (Sum_Sqr_Temp / (n_walkers*n_loops) - (average_temperature * average_temperature))/(n_walkers * n_loops);
        VarToAvrg = (variance/average_temperature);
        if (print_advance) { // print advance percentage every 5% of MAX_NUMBER_OF_LOOPS 
            if (n_loops % advance_step == 0) {
                percentage_advance = (int) ((n_loops*100)/MAX_NUMBER_OF_LOOPS);
                average_temperature = (Sum_Temp) / (n_walkers * n_loops);
                printf("percentage of advance: %d%%  - T_avg: %.6f - variance-to-mean ratio: %.6f \r", 
                    percentage_advance, average_temperature, VarToAvrg); 
            }
        } // close if
    } // close while loop

    // Calculate average temperature
 
    double end = omp_get_wtime();
    // print time of execution
    printf("\n time: %f\n",end-start);
    printf("Average Temperature: %.8f\n", average_temperature);
    printf("Standard Deviation: %.8f\n", sqrt(variance));
    printf("variance-to-mean ratio: %.6f\n", VarToAvrg);
    printf("Walkers: %d \n", (n_walkers*n_loops));
    printf("Loops: %d \n", (n_loops));
    printf("parallel walkers: %d \n", (n_walkers));
   return 0;
 }
