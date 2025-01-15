/******************************************************************************/
/* 
 * Heat equation in 2D with random numbers for (pseudo) steady state.
 * Este ejemplo usa un método de “Random Walk” para estimar la temperatura
 * en un punto del dominio. 
 * 
 * Compilación (ejemplo):
 *   gcc -o heat_WR_transient_BNCH3_v2.exe heat_WR_transient_BNCH3_v2.c ranlib.c rnglib.c -fopenmp
 * Ejecución:
 *   ./heat_WR_transient_BNCH3_v2.exe
 *
 * Referencias (Benchmark):
 *  - 4.3.2 Superposition example 2 (CondBook pag. 108):
 *    https://www.eng.auburn.edu/~dmckwski/mech7210/condbook.pdf
 *  - https://www.professores.uff.br/diomarcesarlobao/wp-content/uploads/sites/85/2017/09/condbook.pdf
 *  - Conduction Heat Transfer Notes for MECH 7210, Daniel W. Mackowski
 *    Auburn University
 *  
 * Observación:
 *   Este ejemplo considera el eje Y invertido respecto a la figura 4.10 
 *   (p.111 del CondBook).
 *
 ******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>

#include "ranlib.h"
#include "rnglib.h"

/******************************************************************************
 *                           DEFINICIÓN DE CONSTANTES                          *
 ******************************************************************************/
#define N_WALKERS 0                   // 0 => usa omp_get_max_threads()
#define MAX_VARIANCE_TO_MEAN_RATIO -1 // -1 => sin control de var/mean
#define DEEP_OF_DOMAIN 1.0
#define WIDTH_OF_DOMAIN 1.0
#define NUMBER_OF_POINTS_IN_X 1000
#define NUMBER_OF_POINTS_IN_Y 1000
#define MAX_NUMBER_OF_LOOPS 1000
#define DEEP_CABLE 1.04
#define SEP_CABLES 0.4
#define TYPE_CABLE_ARRANGEMENT 2
#define BACKFILL_HEIGHT_LAYER 0.4
#define VERBOSE 0
#define VERBOSE2 0
#define VERBOSE3 0
#define VERBOSE4 0
#define PRINT_ADVANCE 1
#define RANDOM_LIB 1
#define TMAX 1000000
#define Factor_dt 1.0

/******************************************************************************
 *                          VARIABLES GLOBALES                                *
 ******************************************************************************/

// Punto donde se calcula la temperatura
// Solucion para x=0.6 ; y=0.4 (y=0.6 en las coordenadas del libro) es 0.14;
double x_point = 0.6; 
double y_point = 0.4; 


// Parámetros principales del dominio
double L = DEEP_OF_DOMAIN;
double W = WIDTH_OF_DOMAIN;
int nx = NUMBER_OF_POINTS_IN_X;
int ny = NUMBER_OF_POINTS_IN_Y;

// Variables de discretización espacial y temporal
double hmin;
int    n_steps;
double d_t_0;
double d_t;

// Datos suelo nativo
double rho_soil1 = 1.0;
double Ce_soil1  = 1.0;
double k_soil1   = 1.0;
double F_soil1   = 0.0;

// Condiciones de borde (Dirichlet)
double Ttop    = 0.0;
double Tbottom = 0.0;
double Tleft   = 0.0;
double Tright  = 0.0;
double Tini    = 0.0;   // Temperatura inicial dentro del dominio

// Condiciones de borde (Neumann)
double qtop    = 0.0;
double qbottom = 0.0;
double qleft   = 0.0;
double qright  = 0.0;

// Tipo de fronteras: 1 => Dirichlet, 2 => Neumann
int b_top    = 1; 
int b_bottom = 2; 
int b_left   = 2; 
int b_right  = 1;

/******************************************************************************
 *                            DECLARACIÓN DE FUNCIONES                        *
 ******************************************************************************/
double  rho   (int i, int j);
double  Ce    (int i, int j);
double  K     (int i, int j);
double  F     (int i, int j, double d_t);
double  U     (int i, int j);

int     calculate_probabilities(int i, int j, double d_t);
double  single_walk(int start_i, int start_j);

int     main(void);

/******************************************************************************
 *                            DEFINICIONES DE FUNCIONES                       *
 ******************************************************************************/

// --------------------------------------------------------------------------
// Densidad en (i, j) - se puede personalizar según cables, suelos, etc.
// --------------------------------------------------------------------------
double rho(int i, int j) {
    (void)i; // no usado
    (void)j; // no usado
    return rho_soil1;
}

// --------------------------------------------------------------------------
// Calor específico en (i, j) - personalizar si hay varias capas
// --------------------------------------------------------------------------
double Ce(int i, int j) {
    (void)i; 
    (void)j; 
    return Ce_soil1;
}

// --------------------------------------------------------------------------
// Conductividad térmica en (i, j)
// --------------------------------------------------------------------------
double K(int i, int j) {
    (void)i;
    (void)j;
    return k_soil1;
}

// --------------------------------------------------------------------------
// Fuente de calor "F" en (i, j), escalar con d_t / (rho*C), etc.
// --------------------------------------------------------------------------
double F(int i, int j, double d_t) {
    // En este ejemplo, se retorna:
    //   1.0 * d_t / (Ce(i,j) * rho(i,j))
    // Podrías reemplazarlo por la lógica real de cables
    (void)i;
    (void)j;
    return 1.0 * d_t / (Ce(i, j) * rho(i, j));
}

// --------------------------------------------------------------------------
// Temperatura inicial en (i, j)
// --------------------------------------------------------------------------
double U(int i, int j) {
    double x = (double)i * hmin;
    double y = (double)j * hmin;
    (void)x;  // si quieres usar x,y en una condición inicial no constante
    (void)y;
    return Tini; 
}

// --------------------------------------------------------------------------
// Calcula la dirección probable del "walker" según la difusividad local
// Retorna un int: 0 => no move, 1 => derecha, 2 => izquierda, 3 => arriba, 4 => abajo
// --------------------------------------------------------------------------
int calculate_probabilities(int i, int j, double d_t) 
{
    // Cantidad de pasos "steps" para la vecindad
    int steps = 1;
    double omega = d_t / (hmin*hmin*steps*steps);

    double K_P = K(i, j);
    double K_W = K(i - steps, j);
    double K_E = K(i + steps, j);
    double K_S = K(i, j + steps);
    double K_N = K(i, j - steps);

    double rho_P = rho(i, j);
    double Ce_P  = Ce(i, j);

    // Probabilidades/cdf 
    double p_P, p_E, p_W, p_N, p_S;
    int verbose = VERBOSE3;

    // Caso: K uniforme
    if (K_P == K_E && K_P == K_W && K_P == K_N && K_P == K_S) {
        double F_c = (K_P / (rho_P * Ce_P)) * omega;
        p_P = 1.0 - 4.0 * F_c; 
        p_E = 1.0 - (3.0 * F_c);
        p_W = 1.0 - (2.0 * F_c);
        p_N = 1.0 - F_c;
        p_S = 1.0; 
    }
    // Caso: K variable
    else {
        float alpha_E = 0.5f * (K_E + K_P) / (rho_P*Ce_P);
        float alpha_W = 0.5f * (K_W + K_P) / (rho_P*Ce_P);
        float alpha_N = 0.5f * (K_N + K_P) / (rho_P*Ce_P);
        float alpha_S = 0.5f * (K_S + K_P) / (rho_P*Ce_P);

        double F_e = alpha_E * omega;
        double F_w = alpha_W * omega;
        double F_n = alpha_N * omega;
        double F_s = alpha_S * omega;
        double F_p = (alpha_W + alpha_E + alpha_N + alpha_S) * omega;

        p_P = 1 - F_p;
        p_E = 1 - (F_p - F_e);
        p_W = 1 - (F_p - F_e - F_w);
        p_N = 1 - (F_p - F_e - F_w - F_n);
        p_S = 1.00;
    }

    // Imprimir debug
    if (verbose) {
        fprintf(stderr,
                "%d: Probabilities: p_P=%.5f p_E=%.5f p_W=%.5f p_N=%.5f p_S=%.5f\n",
                omp_get_thread_num(), p_P, p_E, p_W, p_N, p_S);
    }

    // Generar número aleatorio
    double r;
    if (RANDOM_LIB == 1) {
        r = (double)rand() / (double)RAND_MAX;
    } else {
        // ranlib library
        r = genunf(0.0, 1.0);
    }

    // Seleccionar dirección
    if (r <= p_P) return 0; // no move
    if (r <= p_E) return 1; // right
    if (r <= p_W) return 2; // left
    if (r <= p_N) return 3; // top
    return 4;              // bottom
}

// --------------------------------------------------------------------------
// Simulación de un walker a partir de (start_i, start_j)
// --------------------------------------------------------------------------
double single_walk(int start_i, int start_j)
{
    int i = start_i;
    int j = start_j;
    double temp_walker = 0.0;
    int no_border = 1; 
    int walk_steps = 0;
    double time_ = 0.0;
    double tmax = TMAX;
    int verbose = VERBOSE2;

    d_t = d_t_0; 

    while (no_border == 1 && time_ <= tmax) {
        // Aporte de la fuente de calor
        temp_walker += F(i, j, d_t);
        time_ += d_t;
        walk_steps++;

        int direction = calculate_probabilities(i, j, d_t);

        switch (direction) {
        case 0: 
            // no move
            break;

        case 1: // right
            i++;
            if (i > nx - 1) {
                if (verbose) {
                    fprintf(stderr, "%d: Right Border touched!\n", omp_get_thread_num());
                }
                if (b_right == 1) {
                    // Dirichlet
                    temp_walker += Tright;
                    no_border = 0;
                } else {
                    // Neumann / reflect
                    i = nx - 1;
                }
            }
            break;

        case 2: // left
            i--;
            if (i < 1) {
                if (verbose) {
                    fprintf(stderr, "%d: Left Border touched!\n", omp_get_thread_num());
                }
                if (b_left == 1) {
                    temp_walker += Tleft;
                    no_border = 0;
                } else {
                    i = 1; 
                }
            }
            break;

        case 3: // top
            j--;
            if (j < 1) {
                if (verbose) {
                    fprintf(stderr, "%d: Top Border touched!\n", omp_get_thread_num());
                }
                if (b_top == 1) {
                    temp_walker += Ttop;
                    no_border = 0;
                } else {
                    j = 1;
                }
            }
            break;

        case 4: // bottom
            j++;
            if (j > ny - 1) {
                if (verbose) {
                    fprintf(stderr, "%d: Bottom Border touched!\n", omp_get_thread_num());
                }
                if (b_bottom == 1) {
                    temp_walker += Tbottom;
                    no_border = 0;
                } else {
                    j = ny - 1;
                }
            }
            break;
        }

    } // while (no_border)

    // Si el walker terminó dentro del dominio
    if (no_border == 1) {
        temp_walker += U(i, j);
    }

    if (VERBOSE4) {
        fprintf(stderr,
          "%d: Time=%.2f, i=%d, j=%d, Temp=%.2f, Steps=%d\n",
           omp_get_thread_num(), time_, i, j, temp_walker, walk_steps);
    }

    return temp_walker;
}

// --------------------------------------------------------------------------
// main()
// --------------------------------------------------------------------------
int main(void)
{
    double start_time = omp_get_wtime();

    // Inicializar semilla rand()
    srand(time(NULL));

    // Inicializar ranlib / rnglib
    initialize();
    char phrase[] = "randomizer";
    char buf[64];
    sprintf(buf, "%ld%s", (long)time(NULL), phrase);
    int seed1, seed2;
    phrtsd(buf, &seed1, &seed2);
    set_initial_seed(seed1, seed2);

    // Cálculo de paso hmin
    hmin = L / ((double)nx - 1.0);

    // Cálculo de paso de tiempo d_t_0 (ejemplo con alpha = k/(rho*C))
    double alpha = k_soil1 / (rho_soil1 * Ce_soil1);
    d_t_0 = Factor_dt * (hmin * hmin) / (4.0 * alpha);

    // Punto de inicio del walker
    int start_i = (int)(x_point / hmin);
    int start_j = (int)(y_point / hmin);

    // Configurar número de threads
    int n_walkers = (N_WALKERS == 0) ? omp_get_max_threads() : N_WALKERS;
    omp_set_num_threads(n_walkers);

    printf("Calculation of thermal heat of electric conductors\n");
    printf("Number of parallel walkers: %d\n", n_walkers);
    printf("Max loops: %d\n", MAX_NUMBER_OF_LOOPS);
    printf("Start point: (%d, %d)\n", start_i, start_j);
    printf("hmin: %.6f, L: %.6f\n", hmin, L);

    int n_loops = 0;
    double T_walkers = 0.0;
    double sum_temperature = 0.0;
    double sum_sq_temp = 0.0;
    double Sum_Temp = 0.0;
    double Sum_Sqr_Temp = 0.0;

    double average_temperature = 0.0;
    double variance = 0.0;
    double VarToAvrg = 99999.9;

    // Control del avance
    double perc_advance = 0.0;
    double perc_advance_step = 1.0;
    int    advance_step = (int)(MAX_NUMBER_OF_LOOPS * (perc_advance_step / 100.0));
    int    print_advance = (PRINT_ADVANCE && (advance_step > 0));

    if (print_advance) {
        printf("percentage of advance: %.2f%%  - T_avg: %.6f - var/mean: %.6f \n",
               perc_advance, average_temperature, VarToAvrg);
    }

    // Bucle de iteraciones
    while ((n_loops < MAX_NUMBER_OF_LOOPS) &&
           ((VarToAvrg > MAX_VARIANCE_TO_MEAN_RATIO) || (n_loops < 4))) 
    {
        n_loops++;

        sum_temperature = 0.0;
        sum_sq_temp     = 0.0;

        #pragma omp parallel reduction(+:sum_temperature, sum_sq_temp) private(T_walkers)
        {
            #pragma omp for
            for (int k = 0; k < n_walkers; k++) {
                T_walkers = single_walk(start_i, start_j);
                sum_temperature += T_walkers;
                sum_sq_temp     += (T_walkers * T_walkers);
            }
        } // fin omp parallel

        Sum_Temp     += sum_temperature;
        Sum_Sqr_Temp += sum_sq_temp;

        // Cálculo de media y var
        average_temperature = Sum_Temp / ((double)n_walkers * (double)n_loops);
        double mean_sq  = Sum_Sqr_Temp / ((double)n_walkers * (double)n_loops);
        variance        = mean_sq - (average_temperature * average_temperature);
        variance       /= (n_walkers * n_loops);  // ajusta si tu formula lo requiere

        VarToAvrg = (variance / average_temperature);

        if (print_advance && (n_loops % advance_step == 0)) {
            perc_advance = (double)(n_loops*100.0)/(double)(MAX_NUMBER_OF_LOOPS);
            printf("percentage of advance: %.2f%%  - T_avg: %.6f - var/mean: %.6f\r",
                   perc_advance, average_temperature, VarToAvrg);
            fflush(stdout);
        }
    }

    // Tiempo final
    double end_time = omp_get_wtime();
    double elapsed  = end_time - start_time;

    // Impresión final
    printf("\nDone. Elapsed time: %f s\n", elapsed);
    printf("Average Temperature = %.8f\n", average_temperature);
    printf("Std. Dev = %.8f\n", sqrt(variance));
    printf("Variance-to-Mean ratio = %.6f\n", VarToAvrg);
    printf("Walkers processed = %d\n", (n_walkers * n_loops));
    printf("Loops total = %d\n", n_loops);
    printf("parallel walkers = %d\n", n_walkers);
    printf("delta t_0 = %.8f\n", d_t_0);

    // Guardar resultados en CSV
    FILE *fp = fopen("results_ex3_transient.csv", "a");
    if (fp) {
        fprintf(fp, "%d, %d, %d, %.8f, %.8f, %d, %.8f, %.8f, %d,\n",
                TMAX, start_i, start_j, average_temperature, sqrt(variance),
                (n_walkers*n_loops), hmin, d_t_0, nx);
        fclose(fp);
    }

    return 0;
}
