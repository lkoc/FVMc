// variable refined grid to solve the Heat equation in 2d with random numbers for steady state
//  gcc -o heat_random_v1.exe heat_random_v1.c 

// Main parameters
// lkoc@lkoc-ubnt:~/C_heat$ gcc -o rw_grid rw_grid.c -lm
// lkoc@lkoc-ubnt:~/C_heat$ ./rw_grid


#define DEEP_OF_DOMAIN 10.0 // deep of domain in meters
#define WIDTH_OF_DOMAIN 10.0 // width of domain in meters
#define NUMBER_OF_LEVELS 3 // number of levels of refinement
#define MINIMUM_DISTANCE 0.1 // minimum distance between points in the most refined zone
#define PRINT_GRAPH 1 // 1 to print the graph and 0 to not print the graph

// libraries
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>

int domain2d ( double L, double W, int levels, double hmin) {
// domain2d: creates a 2d domain with a refined grid
// L: deep of domain in meters
// W: width of domain in meters
// levels: number of levels of refinement
// deltaX_far: distance between points in the farthest level of refinement
// h: minimum distance between points in the most refined zone

// for loop to numbers of levels of refinement to create the dominium

int n_levels = 3; 

double h_kX[2*levels -1];
double h_kY[2*levels -1];
int Nx[2*levels-1];
int Ny[2*levels-1];
double Xini[levels+1];
double Xfin[levels];
double Yini[levels+1];
double Yfin[levels];
int Nx_tot = 0;
int Ny_tot = 0;

//Setting of refinement
// X
Xini[0] = 0;
Xini[1] = 0.25*L;
Xini[2] = 0.375*L;
Xini[3] = 0.625*L;
Xfin[0] = Xini[3];
Xfin[1] = 0.75*L;
Xfin[2] = L; //10

// Y
Yini[0] = 0;
Yini[1] = 0.025*W;
Yini[2] = 0.0750*W;
Yini[3] = 0.15*W;
Yfin[0] = Yini[3];
Yfin[1] = 0.75*W;
Yfin[2] = W;

// Levels of refinement
for(int k = 0 ; k <levels; k++) {
    h_kX[k] = hmin * pow(2, levels-1-k); // distance between points
    h_kY[k] = hmin * pow(2, levels-1-k); // distance between points
    Nx[k] = (Xini[k+1] - Xini[k])/h_kX[k]; // number of points
    Ny[k] = (Yini[k+1] - Yini[k])/h_kY[k]; // number of points
    Nx_tot += Nx[k]; // total number of points
    Ny_tot += Ny[k]; // total number of points
}

for(int k = 0 ; k <levels-1; k++) {
    h_kX[k + levels] = hmin * pow(2, k+1); // distance between points;
    h_kY[k + levels] = hmin * pow(2, k+1); // distance between points;;
    Nx[k+levels] = (Xfin[k+1] - Xfin[k])/h_kX[k + levels]; // number of points
    Ny[k+levels] = (Yfin[k+1] - Yfin[k])/h_kY[k + levels]; // number of points
    Nx_tot += Nx[k + levels]; // total number of points
    Ny_tot += Ny[k + levels]; // total number of points
}

// create the 1d vectors with coordinates of x and y per refinement levels and concatenate them
double x[Nx_tot+1];
double y[Ny_tot+1];
int indexX = 0;
int indexY = 0;
x[indexX] = 0;
y[indexY] = 0;
for(int k = 0 ; k < 2*levels -1 ; k++) {
    for (int i=0; i<Nx[k]; i++) {
        indexX++;
        x[indexX] = h_kX[k] + x[indexX-1];
    }
    for (int j=0; j<Ny[k]; j++) {
        indexY++;
        y[indexY] = h_kY[k] + y[indexY-1];
    }
}

// Create a graph showing the 2d grid in png and save the picture in a file called grid.png using gnuplot
if (PRINT_GRAPH) {
    FILE *pipe = popen("gnuplot -persist", "w");
    fprintf(pipe, "set term png\n");
    fprintf(pipe, "set output 'grid.png'\n");
    fprintf(pipe, "set xrange [0:%f]\n", L);
    fprintf(pipe, "set yrange [%f:0]\n", W);
    fprintf(pipe, "set size ratio -1\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set title '2d grid'\n");
    fprintf(pipe, "plot '-' with points pointtype 7 pointsize 0.5 \n");
    for (int i=0; i<Nx_tot; i++) {
        for (int j=0; j<Ny_tot; j++) {
            fprintf(pipe, "%f %f\n", x[i], y[j]);
        }
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    pclose(pipe);
}
// save a vector "x" to a file called gridX.dat with format "i" "x"
// name of variable in the first row and number of points of x in the second row 
FILE *fp;
fp = fopen("gridX.dat", "w");
// print W
fprintf(fp, "%f\n", L);
fprintf(fp, "%d\n", Nx_tot);
fprintf(fp, "i      x\n");
for (int i=0; i<Nx_tot+1; i++) {
    fprintf(fp, "%d %f\n", i, x[i]);
}
fclose(fp);

// save a vector "y" to a file called gridY.dat with format "j" "y"
// name of variable in the first row and number of points of y in the second row 
fp = fopen("gridY.dat", "w");
fprintf(fp, "j      y\n");
fprintf(fp, "%d\n", Ny_tot);
for (int j=0; j<Ny_tot+1; j++) {
    fprintf(fp, "%d %f\n", j, y[j]);
}

return 0;
}

// main function calling the domain2d function with example values assigned to L, W, levels and h 
int main() {
    double L = WIDTH_OF_DOMAIN;
    double W = DEEP_OF_DOMAIN;
    int levels = NUMBER_OF_LEVELS;
    double h = MINIMUM_DISTANCE;
    domain2d(L, W, levels, h);
    return 0;
}

// end of the program