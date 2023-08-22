
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

int main(void)
{
    double start = omp_get_wtime();
    long points = 100000;
    long m = 0;
    unsigned long HAUSNUMERO = 12345;
    double value;
    double x, y;
    char phrase[] = "randomizer";
    char  buf[] = "randomizer1";

    float a;
    float b;
    float high;
    float low;
    int seed1;
    int seed2;

    // Initialize the generators.
    initialize ( );
        // use sprintf() function to add string "phrase" to time and convert to a string in buf
    sprintf ( buf, "%d%s", time(NULL), phrase );

    // Set the seeds based on the phrase.
    phrtsd ( buf, &seed1, &seed2 );
    // Initialize all generators.
    set_initial_seed ( seed1, seed2 );

    // Select the parameters at random within a given range.
    low = 0.0;
    high = 1.0;
    
    int threads = omp_get_max_threads();
    omp_set_num_threads(threads);
    
    #pragma omp parallel reduction (+: m ) private (x, y) 
        {
        #pragma omp for nowait
        for(int i = 0; i < points; i++)
            {
            x = genunf ( low, high );
            y = genunf ( low, high );
            m += ( x * x + y * y <= 1.0);
            }
        }
    double end = omp_get_wtime();
    //print of number of points
    printf("points: %ld\n", points);
    printf("time: %f\n",end-start);
    printf("Pi is roughly %lf\n", (double) 4*m / (double) points);
 }

