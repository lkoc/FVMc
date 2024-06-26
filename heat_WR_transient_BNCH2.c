// Heat equation in 2d with random numbers for steady state
//  gcc -o heat_WR_transient_BNCH2.exe heat_WR_transient_BNCH2.c  ranlib.c rnglib.c -fopenmp
// ./heat_WR_transient_BNCH2.exe
//
//
// Benchmark 1: ejemplo 2 del paper:
// A meshless model for transient heat conduction in functionally graded materials
//
//  https://www.researchgate.net/publication/226961451_A_meshless_model_for_transient_heat_conduction_in_functionally_graded_materials
//
// Wang, H., Qin, Q. H., & Kang, Y. (2006). A meshless model for transient heat 
// conduction in functionally graded materials. Computational mechanics, 38(1), 51-60.

#define PI 3.14159265358979323846

// Main parameters
#define N_WALKERS 0 // number of parallel walkers 0 for maximum hardware number of threads
#define MAX_VARIANCE_TO_MEAN_RATIO -1 //1.0 //0.025 //-1.0 // maximum variance to mean ratio to stop the simulation (-1 for no control by variance-to-mean ratio)
#define DEEP_OF_DOMAIN 0.04 // deep of domain in meters
#define WIDTH_OF_DOMAIN 0.04 // width of domain in meters
#define NUMBER_OF_POINTS_IN_X 200 // 7500
#define NUMBER_OF_POINTS_IN_Y 200 //7500
#define MAX_NUMBER_OF_LOOPS 2000 // maximum number of loops
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
#define TMAX 30.0 // maximum temperature of the domain 
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
// Point to calculate the temperature
 ///////////////////////////////////////
 ///////////////////////////////////////
double x_point = 0.02; //x_cable21;          //   L/2; // x coordinate of the point
double y_point = 0.02; //y_cable21 ; // - cable_diam/2-0.002; // y coordinate of the point
double lamb = 50.0; // parameter of K function
// Solucion para lamb=50, x=0.02 ; y=0.02 es 0.68@ t=10 ; 0.73@t=30

// Main parameters
double L = DEEP_OF_DOMAIN; // length of domain
double W = WIDTH_OF_DOMAIN; // width of domain
int nx= NUMBER_OF_POINTS_IN_X; // number of points in x direction
int ny= NUMBER_OF_POINTS_IN_Y; // number of points in y direction

double hmin; // spatial step
int n_steps; // number of steps in each walk
double d_t_0; // initial time step
double d_t; // variable time step

// Native Soil thermal data for calculations
double rho_soil1 = 1000.0; // density of soil_1 in kg/m3 - fuente: https://www.cableizer.com/documentation/zeta_soil/
double Ce_soil1 = 1000.0; // specific heat of soil 1 in J/(kg K) - fuente: https://www.cableizer.com/documentation/c_p_soil/
double k_soil1 = 17.0;// 0.2; // thermal conductivity of soil 1 in W/(m.K) (resistivity = 5)- fuente: https://link.springer.com/article/10.1007/s10765-022-03119-5
//double F_soil1 = 0.0; // heat source in soil 1 W/m3

//boundaries temperatures and gradients
// Dirichlet boundary conditions
double Ttop    = 0.0; // temperature at the top of the domain
double Tbottom = 0.0; // temperature at the bottom of the domain
double Tleft   = 0.0; // temperature at the left of the domain
double Tright  = 1.0; // temperature at the right of the domain
double Tini    = 0.0; // initial temperature of the domain

// Neumann boundary conditions
double qtop    = 0.0; // heat flux at the top of the domain
double qbottom = 0.0; // heat flux at the bottom of the domain
double qleft   = 0.0; // heat flux at the left of the domain
double qright  = 0.0; // heat flux at the right of the domain

// kind of boundaries
int b_top    = 2; // 1: Dirichlet, 2: Neumann
int b_bottom = 2; // 1: Dirichlet, 2: Neumann
int b_left   = 1; // 1: Dirichlet, 2: Neumann
int b_right  = 1; // 1: Dirichlet, 2: Neumann


///////////////////////////////////////////////////////////////////////////////
// Function to calculate density at a given coordinate (i, j)               //
///////////////////////////////////////////////////////////////////////////////
double rho(int i, int j) {
    // Define density calculation based on (x, y) coordinates
    return rho_soil1;
}

///////////////////////////////////////////////////////////////////////////////
// Function to calculate specific heat Ce at a given coordinate (i, j)      //
///////////////////////////////////////////////////////////////////////////////
double Ce(int i, int j) {
    // Define Specific Heat calculation based on (x, y) coordinates
    return Ce_soil1;
}

///////////////////////////////////////////////////////////////////////////////
// Function to calculate thermal conductivity "K" at a coordinate (i, j)     //
///////////////////////////////////////////////////////////////////////////////
double K(int i, int j) {
    // Define thermal conductivity calculation based on (x, y) coordinates
    return k_soil1 * exp(lamb * i * hmin);
}

///////////////////////////////////////////////////////////////////////////////
// Function to calculate heat source "F" at a given coordinate (i, j)          //
///////////////////////////////////////////////////////////////////////////////
double F(int i, int j, double d_t) {
    // Define source heat calculation based on (x, y) coordinates
    return 0.0*d_t/(Ce(i,j) * rho(i,j)); 
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
int calculate_probabilities( int i, int j,  double d_t) {
    // Define the probabilities of the walker to move to the next position
    // based on the thermal conductivity of the material at the next position
    double Ce_E, Ce_W, Ce_S, Ce_N , Ce_P;    // Coefficients of the temperature terms of the heat equation at interior node
    double p_P, p_E, p_W, p_N, p_S;      // Probabilities of the walker to move to the next position
    double K_E, K_W, K_N, K_S, K_P; // thermal conductivity around the point P --> (i,j)
    double F_e, F_w, F_n, F_s, F_p; // linearly interpolated values of thermal conductivity
    double F_c ; // linearly interpolated values of thermal conductivity of central point    
    double rho_E, rho_W, rho_N, rho_S, rho_P; // linearly interpolated values of thermal conductivity
    int steps =1;
    double omega = d_t / (1.0*hmin*hmin*steps*steps); //1.0*d_t / (hmin*hmin*steps*steps); // time step divided by h^2
    
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
        // Valores de alpha
        float alpha_E = 0.5*(K_E + K_P) / (rho_P*Ce_P);
        float alpha_W = 0.5*(K_W + K_P) / (rho_P*Ce_P);
        float alpha_N = 0.5*(K_N + K_P) / (rho_P*Ce_P);
        float alpha_S = 0.5*(K_S + K_P) / (rho_P*Ce_P);

        // Valores de F (coeficientes) interpolados linealmente
        
        F_e = alpha_E * omega ; // valor de F entre E y P
        F_w = alpha_W * omega ; // valor de F entre W y P
        F_n = alpha_N * omega ; // valor de F entre N y P
        F_s = alpha_S * omega ; // valor de F entre S y P
        F_p = (alpha_W + alpha_E + alpha_N + alpha_S) * omega ; // valor de F en P

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

    double time = 0.0; // time point for the walker
    double tmax = TMAX; // maximum time for the walker
    
    d_t = d_t_0; // time step size for cooper conductivity
    // loop over the walk length 
    while (no_border == 1 && time <= tmax) {
        temp_walker += F(i,j,  d_t)  ; // *steps;
        time += d_t; // update the time
        
        walk_steps++; // update the walk steps number
        //Add heat source to the temperature
        
        //Calculate the probabilities and move the walker 
        direction = calculate_probabilities(i, j, d_t); // retorna direccion random y actuliza paso de tiempo
        switch (direction) {
            case 0:
                //no_border = 0; // No border
                //temp_walker += U(i,j);
                break; // No move
            case 1:
                i = i + 1; // move to the right
                if (i > nx - 1) { // Right border
                    if (verbose == 1) {fprintf(stderr," %d: Right Border touched! \n ", omp_get_thread_num());}
                    if (b_right == 1) { // Dirichlet border
                        temp_walker += Tright; 
                        no_border = 0; // No border
                    } else { // Neumann border (with q=0)
                        i= nx-1; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 2:
                i = i - 1; // move to the left
                if (i < 1) {  // Left border
                    if (verbose == 1) {fprintf(stderr," %d: Left Border touched! \n ", omp_get_thread_num());}
                    if (b_left == 1) { // Dirichlet border
                        temp_walker += Tleft;
                        no_border = 0; // No border 
                    } else { // Neumann border  (with q=0)
                        i =1; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    }
                }
                break;
            case 3:
                j = j - 1; // move to the top
                if (j < 1) { // Top border
                    if (verbose == 1) {fprintf(stderr," %d: Top Border touched! \n ", omp_get_thread_num());}
                    if (b_top == 1) { // Dirichlet border
                        temp_walker += Ttop;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= 1; // reflective case: move the walker inside the domain
                        no_border = 1; // inside the domain
                    } 
                }
                break;
            case 4:
                j = j + 1; // move to the bottom
                if (j > ny - 1) { // Bottom border
                    if (verbose == 1) {fprintf(stderr," %d: Bottom Border touched! \n ", omp_get_thread_num());}
                    if (b_bottom == 1) { // Dirichlet border
                        temp_walker += Tbottom;
                        no_border = 0; // No border
                    } else { // Neumann border  (with q=0)
                        j= ny-1; // reflective case: move the walker inside the domain
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

// REsultado teórico del paper mencionado en el encabezado BNCmark 1
double u_teor(double t, double xvec) {
    double u = 0.0;
    for (int i = 0; i < 10000; i++) {
        double mu = (2.0 * i + 1.0) * PI / 2.0;
        u += pow(-1.0, i) * 4.0 / ((2.0 * i + 1.0) * PI) *
              cos(mu * xvec) * exp(-pow(mu, 2.0) * t);
    }
    return 1.0 - u;
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
    double alpha = K(nx, 0) / (rho_soil1 * Ce_soil1);
    d_t_0 = Factor_dt * (hmin*hmin )  / (4 * alpha); // time step size for cooper conductors
    int t_max = TMAX ;

    ///////////////////////////////////////
    ///////////////////////////////////////
    // Initial point of walkers
    int start_i = (int) (x_point/hmin); // starting x coordinate
    int start_j = (int) (y_point/hmin); // starting y coordinate
    
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
    
    printf("Number of parallel walkers: %d\n", n_walkers);
    printf("Number of loops: %d\n", MAX_NUMBER_OF_LOOPS);
    // print starting point of walkers
    printf("Starting point of walkers: (%d,%d)\n", start_i, start_j);
    // print h and L
    printf("hmin: %.6f\n", hmin);
    
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
    printf("Tem teor: %f\n", u_teor(TMAX, x_point));

    // write in a single text row Average Temperature, standard deviation and  n_walkers*n_loops to file to a file (new or add to existent) 
    FILE *fp; // file pointer
    fp = fopen("results_ex3_transientBNCH2.csv", "a"); // open file to append
    fprintf(fp, "%d, %d, %d, %.8f, %.8f, %d, %.8f, %.8f, %d, \n", t_max, start_i, start_j, average_temperature, sqrt(variance), n_walkers*n_loops, hmin, d_t_0 ,nx) ;
    fclose(fp);
   return 0;
}
