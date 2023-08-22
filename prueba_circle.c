// libraries
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
///////////////////////////////////////////////////////////////////////////////
// Function to check if coordinates (x,y) are inside a circle of radio r     //
// and center (x0,y0)                                                       //
///////////////////////////////////////////////////////////////////////////////
int circle(double x, double y, double x0, double y0, double r) {
    return (((x - x0) * (x - x0) + (y - y0) * (y - y0)) <= (r * r));
    }

///////////////////////////////////////////////////////////////////////////////
// main function  with two external arguments x and y                        //
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

// convert external arguments *argv[] to double type variables x and y 
double x=-1;
double y=-1;

if (argc >= 2) {
    x = atof(argv[1]);
    y = atof(argv[2]);
    }

double r=0.05; //radio del circulo
double x0=1.0; //centro del circulo
double y0=1.0; //centro del circulo

//call function circle and prit if the (x,y) is inside or outside the circle
if (circle(x,y,x0,y0,r)==1){
    printf("The point (%f,%f) is inside the circle\n",x,y);
    }
else{
    printf("The point (%f,%f) is outside the circle\n",x,y);
    }
return 0;
}




