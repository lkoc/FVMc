// Thermal calculation for a three circuits of power cables in trefoil in direct underground buried
// random walk simulation of the temperature in the variable soil and in the cable
// the soil is considered non-homogeneous and non-isotropic
// the cable is considered homogeneous and isotropic
// the soil is considered as a 2D domain
// 
// The dinamic random walk finish when the time os ovre in the time sptep is variable in order to maximace the value of probabiliries and
// acelerate tha sampling of the random walk 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <omp.h>
//#include "mt19937ar.h"

// Main parameters
#define N_WALKERS 4 // number of parallel walkers
#define MAX_VARIANCE_TO_MEAN_RATIO 0.1 // maximum variance to mean ratio to stop the simulation
#define DEEP_OF_DOMAIN 2.0 // deep of domain in meters
#define WIDTH_OF_DOMAIN 2.0 // width of domain in meters
#define NUMBER_OF_POINTS_IN_X 20000
#define NUMBER_OF_POINTS_IN_Y 20000
#define MAX_TIME 1000000 // maximum time in seconds
#define MAX_NUMBER_OF_LOOPS 1000000 // maximum number of loops
#define DEEP_CABLE 1.0  // deep of cable in meters
#define SEP_CABLES 0.2  // separation between cables in meters between trifoils
#define VERBOSE 0 // 1 to print detailed results in the screen
#define VERBOSE2 0 // 1 to print loops results in the screen
#define VERBOSE3 0 // 1 to print when heat source is activated
#define TIME_SCALE 0.99// time scale to accelerate the simulation

// max and min functions
//////////////////////////////////////////////////////////////////////////////////
#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})
//////////////////////////////////////////////////////////////////////////////////

// Global Variables

// Main parameters
double L = DEEP_OF_DOMAIN; // length of domain
double W = WIDTH_OF_DOMAIN; // width of domain
int nx= NUMBER_OF_POINTS_IN_X; // number of points in x direction
int ny= NUMBER_OF_POINTS_IN_Y; // number of points in y direction
double Tmax= MAX_TIME; // maximum time

double h; // spatial step
double delta_t; // temporal step

int n_walkers = N_WALKERS; // number of walkers
int n_steps; // number of steps in each walk

// Aluminum 630 mm2 data for calculations
double cable_area= 630 ;// cross area of cable 630 mm2
double cable_diam = 0.0552; // diameter of cable 630 mm2 - fuente: http://es.lemeicable.com/product/108.html
double cond_diam =  0.0301; // diameter of conductor 630 mm2 - fuente: http://es.lemeicable.com/product/108.html 
double res_cond = 0.0000601; // ohm/m  - electric resistivity of cable 630 mm2 0.0469 ohm/km at 20°C DC - fuente: http://es.lemeicable.com/product/108.html
double rho_cond = 2700.0; // density of aluminum cable in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_cond = 900.0; // specific heat of aluminum cable in J/(kg K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double k_cond = 237.0; // thermal conductivity of aluminum cable in W/(m.K) - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html

// datos XLPE
double rho_xlpe = 920.0; // density of XLPE in kg/m3 - fuente: http://www.engineeringtoolbox.com/metal-alloys-densities-d_50.html
double Ce_xlpe = 2174.0; // specific heat of XLPE in J/(kg K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
double k_xlpe = 0.28; // thermal conductivity of XLPE in W/(m.K) - fuente: https://www.emworks.com/application/numerical-modeling-of-the-thermal-behavior-of-xlpe-power-cable
// double F_xlpe = 0.0; // heat source in cable W/m3

// Native Soil thermal data for calculations
double rho_soil1 = 2000.0; // density of soil_1 in kg/m3 - fuente: https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil1 = 800.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil1 = 0.125; // thermal conductivity of soil 1 in W/(m.K) (resistivity = 8)- fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5
// double F_soil1 = 0.0; // heat source in soil 1 W/m3

// Backfill Soil thermal data for calculations
double rho_soil2 = 1500.0; // density of soil_1 in kg/m3 - fuente:  https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil2 = 900.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil2 = 0.33; // thermal resistivity of soil 1 in W/(m.K) (resistivity = 3) - fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5 
double soil2_layer = 0.25; // height of backfill soil layer in m
// double F_soil2 = 0.0; // heat source in soil 1 W/m3

// max, min values
double rho_min = 1500; // minimum density
double Ce_min = 900; // minimum specific heat
double k_max = 237; // maximum thermal diffusivity

//Location of conductors (x=0 is the left side of the domain, y=0 is the top of the domain)
double deep_cable = DEEP_CABLE; // depth of cable in m, from the surface to the cable touch the soil
double separation_cable = SEP_CABLES; // separation between cable circuits in m

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
double Tbottom = 20.0; // temperature at the bottom of the domain
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

// Function to check if coordinates (x,y) are inside a circle of radio r and center (x0,y0)
int circle(double x, double y, double x0, double y0, double r) {
    return (((x - x0) * (x - x0) + (y - y0) * (y - y0)) <= (r * r));
    }

// Function to calculate density at a given coordinate (i, j)
double rho(int i, int j) {
    // Define density calculation based on (x, y) coordinates
    double rho_ = rho_soil1;
    
    double x = i*h;
    double y = j*h;

    // (x,y) inside soil 2 area - Backfill Soil
    if (y < (deep_cable - cable_diam*(0.5774 + 0.5) + soil2_layer/2) &&
        y > (deep_cable - cable_diam*(0.5774 + 0.5) - soil2_layer/2)) {
        rho_ = rho_soil2;
   }

    // Circuit 1, cable 11 , cable 12 and cable 13   - Central circuit
    if (circle(x, y, x_cable11, y_cable11, cable_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cable_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cable_diam/2)) {
        rho_ = rho_xlpe;
    }
    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        rho_ = rho_cond;
    }

    // Circuit 2, cable 21 , cable 22 and cable 23 - Left circuit
    if (circle(x, y, x_cable21, y_cable21, cable_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cable_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cable_diam/2)) {
        rho_ = rho_xlpe;
    }
    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        rho_ = rho_cond;
    }
    
    // Circuit 3, cable 31 , cable 32 and cable 33 - Right circuit
    if (circle(x, y, x_cable31, y_cable31, cable_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cable_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cable_diam/2)) {
        rho_ = rho_xlpe;
    }

    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        rho_ = rho_cond;
    }
    return rho_;
}

// Function to calculate specific heat Ce at a given coordinate (i, j)
double Ce(int i, int j) {
    // Define density calculation based on (x, y) coordinates
    double Ce_ = Ce_soil1;
    
    double x = i*h;
    double y = j*h;

    // (x,y) inside soil 2 area - Backfill Soil
    if (y < (deep_cable - cable_diam*(0.5774 + 0.5) + soil2_layer/2) &&
        y > (deep_cable - cable_diam*(0.5774 + 0.5) - soil2_layer/2)) {
        Ce_ = Ce_soil2;
   }

    // Circuit 1, cable 11 , cable 12 and cable 13   - Central circuit
    if (circle(x, y, x_cable11, y_cable11, cable_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cable_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cable_diam/2)) {
        Ce_ = Ce_xlpe;
    }
    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        Ce_ = Ce_cond;
    }

    // Circuit 2, cable 21 , cable 22 and cable 23 - Left circuit
    if (circle(x, y, x_cable21, y_cable21, cable_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cable_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cable_diam/2)) {
        Ce_ = Ce_xlpe;
    }
    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        Ce_ = Ce_cond;
    }
    
    // Circuit 3, cable 31 , cable 32 and cable 33 - Right circuit
    if (circle(x, y, x_cable31, y_cable31, cable_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cable_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cable_diam/2)) {
        Ce_ = Ce_xlpe;
    }

    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        Ce_ = Ce_cond;
    }
    return Ce_;
}

// Function to calculate thermal conductivity "k" at a given coordinate (i, j)
double k(int i, int j) {
    // Define thermal conductivity calculation based on (x, y) coordinates
    double k_ = k_soil1;
    
    double x = i*h;
    double y = j*h;

    // (x,y) inside soil 2 area - Backfill Soil
    if (y < (deep_cable - cable_diam*(0.5774 + 0.5) + soil2_layer/2) &&
        y > (deep_cable - cable_diam*(0.5774 + 0.5) - soil2_layer/2)) {
        k_ = k_soil2;
   }

    // Circuit 1, cable 11 , cable 12 and cable 13   - Central circuit
    if (circle(x, y, x_cable11, y_cable11, cable_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cable_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cable_diam/2)) {
        k_ = k_xlpe;
    }
    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        k_ = k_cond;
    }

    // Circuit 2, cable 21 , cable 22 and cable 23 - Left circuit
    if (circle(x, y, x_cable21, y_cable21, cable_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cable_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cable_diam/2)) {
        k_ = k_xlpe;
    }
    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        k_ = k_cond;
    }
    
    // Circuit 3, cable 31 , cable 32 and cable 33 - Right circuit
    if (circle(x, y, x_cable31, y_cable31, cable_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cable_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cable_diam/2)) {
        k_ = k_xlpe;
    }

    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        k_ = k_cond;
    }
    return k_;
}

// Function to calculate heat source at a given coordinate (i, j)
double F(int i, int j) {
    // Define source heat calculation based on (x, y) coordinates
    double F_ = 0.0; // Assumed that only conductors generate heat
    
    double x = i*h;
    double y = j*h;

    if (circle(x, y, x_cable11, y_cable11, cond_diam/2) || 
        circle(x, y, x_cable12, y_cable12, cond_diam/2) || 
        circle(x, y, x_cable13, y_cable13, cond_diam/2)) {
        F_ = F_c1;
    }

    if (circle(x, y, x_cable21, y_cable21, cond_diam/2) || 
        circle(x, y, x_cable22, y_cable22, cond_diam/2) || 
        circle(x, y, x_cable23, y_cable23, cond_diam/2)) {
        F_ = F_c2;
    }
    
    if (circle(x, y, x_cable31, y_cable31, cond_diam/2) || 
        circle(x, y, x_cable32, y_cable32, cond_diam/2) || 
        circle(x, y, x_cable33, y_cable33, cond_diam/2)) {
        F_ = F_c3;
    }
    return F_;
}

// Function to return initial temperature given coordinate (x, y)
double U(int i, int j) {
    double x = i*h;
    double y = j*h;
    // Define heat source calculation based on (x, y) coordinates
    return Tini; // Initial temperature
}  

// Border function must to return two variables: the border value and a boolean indicating if the walker is out of bounds
int border_temp(int* i, int* j, double* temp_sum) {  // return 0=no border, 1=border
    // Define different border temperatures here
    int verbose = VERBOSE ; // 1 = print border touched
    if (*i < 1) {  // Left border
        if (verbose == 1) {
            printf("Border touched! \n ");
            // print walker position
            printf("i = %d, j = %d\n", *i, *j);
            //exit(2);
        }
        if (b_left == 1) { // Dirichlet border
            *temp_sum += Tleft;
            return 0; // No border
         } else { // Neumann border
            *i++; // reflective case: move the walker inside the domain
            return 1; // inside the domain
            }
        } 
        else if (*i > nx - 1) { // Right border
            if (verbose == 1) {printf("Border touched!");}
            if (b_right == 1) { // Dirichlet border
                *temp_sum += Tright;
                return 0; // No border
            } else { // Neumann border
                *i--; // reflective case: move the walker inside the domain
                return 1; // inside the domain
            } 
        }
        else if (*j < 1) { // Top border
            if (verbose == 1) {printf("Border touched!");}
            if (b_top == 1) { // Dirichlet border
               *temp_sum += Ttop;
               return 0; // No border
            } else { // Neumann border
                *j++; // reflective case: move the walker inside the domain
                return 1; // inside the domain
            } 
        }
        else if (*j > ny - 1) { // Bottom border
            if (verbose == 1) {printf("Border touched!");}    
            if (b_bottom == 1) { // Dirichlet border
                *temp_sum += Tbottom;
                return 0; // No border
            } else { // Neumann border
                *j--; // move the walker inside the domain
                return 1; // inside the domain
            }
        }

    return 1; // default case (inside the domain)
}

int calculate_probabilities( int i, int j) {
    // calculate probabilities based on position and calculation parameters
    // x, y
    int verbose = VERBOSE; // 1 = print probabilities, 0 = don't print probabilities
    delta_t = h*h * Ce(i,j) *rho(i,j)/(4*k(i,j)); // maximum time step
    
    double Kxy =   k(i     , j);
    double Kxm_1y = k(i - 1 , j);
    double Kxp_1y = k(i + 1 , j);
    double Kxym_1 = k(i     , j - 1 );
    double Kxyp_1 = k(i     , j + 1);
    
    double Cxyrhoxy = Ce(i, j) * rho(i, j);

    double halfKxym1Kxy = 0.5 * (Kxym_1 + Kxy);
    double halfKxyp1Kxy = 0.5 * (Kxyp_1 + Kxy);

    double halfKxm1yKxy = 0.5 * (Kxm_1y + Kxy);
    double halfKxp1yKxy = 0.5 * (Kxp_1y + Kxy);

    // Variable time step to ensure stability and reduce sampling 
    delta_t= TIME_SCALE * h*h * Cxyrhoxy/(4*k(i,j)); // maximum time step
    double s= delta_t;
    // Return the probabilities as an array [Po, Po+Pe, Po+Pe+Pw, Po+Pe+Pw+Pn , Po+Pe+Pw+Pn+Ps]
    double p_o , p_oe, p_oew, p_oewn , p_oewns ;
    p_o   = 1 - 4*s*Kxy/(h*h*Cxyrhoxy);
    p_oe  = 1 - s/(h*h*Cxyrhoxy)*(3*Kxy - halfKxm1yKxy + halfKxp1yKxy);
    p_oew = 1 - 2*s*Kxy/(h*h*Cxyrhoxy);
    p_oewn = 1 - s/(h*h*Cxyrhoxy)*(Kxy - halfKxym1Kxy + halfKxyp1Kxy);
    p_oewns = 1;

    // print p_0, P_1, P_2, P_3, P_4 and the sum of all probabilities
    if (verbose == 1) {
        printf("p_o = %f, p_oe = %f, p_oew = %f, p_oewn = %f, p_oewns = %f \n", p_o, p_oe, p_oew, p_oewn, p_oewns);
    }

    // random number between 0 and 1
    double r =  (double)rand() / (double)RAND_MAX; //   (double)rand() / (double)RAND_MAX;
    if (verbose == 1) {printf(" random number: %.6f\n", r);}
    // random directions or walker
    if (r <= p_o) { // stay in the same position
        return 0;
    }
    if (r <= p_oe) { // move to the right
        return 1;
    }
    if (r <= p_oew) { // move to the left
        return 2;
    }
    if (r <= p_oewn) { // move to the top
        return 3;
    }

    return 4; // move to the bottom
}

double single_walk(int start_i, int start_j) {
    // single random walk
    int i = start_i;
    int j = start_j;

    double ttime = 0;
    int direction;
    double temp_sum = 0 ; //U(i,j);
    int no_border = 1;
    int verbose = VERBOSE ; // 1 = print walker position, 0 = don't print walker position
    int kk=0;

    // loop over the walk length 
    while (no_border == 1 && ttime < Tmax) {
        kk++;
        //Add heat source to the temperature
        temp_sum += F(i,j) * delta_t/(Ce(i,j)*rho(i,j)); 
        if (F(i,j) > 0.0 ) {
            if (VERBOSE3 == 1) {
                printf("Heat source added!\n");
                printf("temp sum: F(i,j): %.8f\n", temp_sum);}}
        if (verbose == 1) {
            printf("i: %d, j: %d\n", i, j);
            printf(" Source: %.2f\n", F(i,j) * delta_t/(Ce(i,j)*rho(i,j)));
            printf(" Time: %.2f\n", ttime);
        }
        
        //Call border_temp function and check if the walker is out of bounds when temp is different from 9999.9
        no_border = border_temp(&i, &j, &temp_sum);

        //Calculate the probabilities and move the walker 
        direction = calculate_probabilities(i, j);
        ttime += delta_t;
        switch (direction) {
            case 0:
                temp_sum += U(i,j); // walker stay in the same position and tally the temperature of previus time step
                no_border= 0; // stop the walk because the walker remain in the same position
                break;
            case 1:
                i++; // move to the right
                break;
            case 2:
                i--; // move to the left
                break;
            case 3:
                j++; // move to the top
                break;
            case 4:
                j--; // move to the bottom
                break;
        }

    }

    if (no_border == 1) {
        temp_sum += U(i,j); 
       }
    if (verbose == 1) {printf("Walker sum Temp: %.2f\n", temp_sum);}
    return temp_sum;
}

int main() {
    // space step size h
    h= L/(nx-1); // nx-1 because we have nx points in the domain and nx-1 intervals

    // Initialize values for circuits
    //Circuit 1 - three phase arrangement in trifoil - Central circuit
    x_cable11 = L/2; // x coordinate of cable 1 , central and upper cable
    y_cable11 = deep_cable - cable_diam*(0.866  + 0.5); // y coordinate of cable 1
    x_cable12 = x_cable11 - cable_diam; // x coordinate of cable 2, left and down cable
    y_cable12 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable13 = x_cable11 + cable_diam; // x coordinate of cable 3, right and down cable
    y_cable13 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // current in circuit 1
    double Ic1 = 250.0; // current in cable
    F_c1 = Ic1*Ic1*res_cond/(cable_area/1000000); // heat source in cable W/m3

    //Circuit 2 - three phase arrangement in trifoil - Left circuit
    x_cable21 = x_cable11 - separation_cable; // x coordinate of cable 1 , central and upper cable
    y_cable21 = y_cable11; // y coordinate of cable 1
    x_cable22 = x_cable21 - cable_diam; // x coordinate of cable 2, left and down cable
    y_cable22 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable23 = x_cable21 + cable_diam; // x coordinate of cable 3, right and down cable
    y_cable23 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // current in circuit 2
    double Ic2 = 250.0; // current in cable
    F_c2 = Ic2*Ic2*res_cond/(cable_area/1000000); // heat source in cable W/m3

    //Circuit 3 - three phase arrangement in trifoil - Right circuit
    x_cable31 = x_cable11 + separation_cable; // x coordinate of cable 1 , central and upper cable
    y_cable31 = y_cable11; // y coordinate of cable 1
    x_cable32 = x_cable31 - cable_diam; // x coordinate of cable 2, left and down cable
    y_cable32 = deep_cable - cable_diam/2; // y coordinate of cable 2
    x_cable33 = x_cable31 + cable_diam; // x coordinate of cable 3, right and down cable
    y_cable33 = deep_cable - cable_diam/2; // y coordinate of cable 3
    
    // current in circuit 3
    double Ic3 = 250.0; // current in cable
    F_c3 = Ic3*Ic3*res_cond/(cable_area/1000000); // heat source in cable W/m3

    // Point to calculate the temperature
    double x_point = x_cable11; //   L/2; // x coordinate of the point
    double y_point = y_cable11; //  deep_cable - 1.15*cable_diam; // y coordinate of the point
    // Initial point of walkers
    int start_i = (int) (x_point/h); // starting x coordinate
    int start_j = (int) (y_point/h); // starting y coordinate

    clock_t start, end;
    // Start clock
    start = clock();
     
//    double smax = h*h * Ce_min*rho_min/(4*k_max); // maximum time step
//    //delta_t = min((2*Tmax/(nx*nx)), smax*0.99);// time step
//    delta_t = smax*0.95; // minimize the p that remain in the same postion
//    Tmax = delta_t * L*L/(4*h*h); // maximum time 
//    n_steps= (int) (Tmax/delta_t+0.5); // number of steps in each walk

    // Seed random number generator
    srand(time(NULL));

    // Variables for calculating average and standard deviation
    double sum_temperature = 0.0;
    double sum_squared_temperature = 0.0;
    double average_temperature = 9999.0;
    double variance = 9999.0; 
    int n_loop = 0; // Number of loops
    
    // Outer loop to control accumulated standard deviation
    //while ((variance /average_temperature > MAX_VARIANCE_TO_MEAN_RATIO ||(n_loop < MAX_NUMBER_OF_LOOPS/10) ) && (n_loop < MAX_NUMBER_OF_LOOPS)) {
    while (n_loop < MAX_NUMBER_OF_LOOPS) {
        // Increment loop counter
        n_loop++; 
        // Parallel loop with reduction for sum and sum of squares
        #pragma omp parallel for reduction(+:sum_temperature) reduction(+:sum_squared_temperature)
        for (int k = 0; k < n_walkers; k++) {
            double T_walker = single_walk(start_i, start_j);
            sum_temperature += T_walker;
            sum_squared_temperature += T_walker*T_walker;
            //printf(" T Walker: %.2f\n", T_walker);
            //printf(" T average: %.2f\n", sum_temperature/(n_loop*(k+1)));
        }
        // Calculate average
        average_temperature = (sum_temperature) / (n_walkers * n_loop);
        
        // Verbose output
        if (VERBOSE2 == 1) {
            printf(" T average: %.6f\n", average_temperature);
        }
        
        // Calculate variance
        variance = (sum_squared_temperature / (n_walkers*n_loop) - (average_temperature * average_temperature))/(n_walkers * n_loop);
    }
    // End clock
    end = clock();
 
    // Calculate the time taken for the execution
    double time_taken = ((double)(end - start))/CLOCKS_PER_SEC;
    
    // Print the Calculated execution time
    printf("Execution time:%.4f seconds \n", time_taken);
    printf("Average Temperature: %.8f\n", average_temperature);
    printf("Standard Deviation: %.8f\n", sqrt(variance));
    printf("variance-to-mean ratio: %.6f\n", variance/average_temperature);
    printf("Walkers: %d \n", (n_walkers*n_loop));
    printf("Loops: %d \n", (n_loop));
   return 0;
}