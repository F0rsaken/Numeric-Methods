#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <math.h>
// #include "../../lib/matrix-lib.h"
// #include "../../lib/interpolation-lib.h"
#include "../../lib/diff-equations-lib.h"
using namespace std;


/**
 * Zadana funkcja:
 * y' - kmy sin(mx) = k^2 m sin(mx) cos(mx)
 * y(x0) = a
 * 
 * Zadane parametry:
 * x0 = pi / 2
 * xk = (7*pi)/ 2
 * m = 1
 * k = 5
 * 
 * Finalna postać:
 * y' - 5y*sin(x) = 25*sin(x)*cos(x)
 * y' = 25*sin(x)*cos(x) + 5y*sin(x) = 5sin(x) * ( 5cos(x) + y )
 * y' = 5sin(x) * (5cos(x) + y)
 * 
 * a = TODO: wyliczyć
 * 
 * rozwiazanie:
 *  y(x) = e^(-5 * cos(x)) - 5*cos(x) + 1
 */

/**
 * Argumenty:
 *  - h
 *  - 1 | 2 (Euler | Runge-Kutta)
 **/

double x0 = M_PI/2;
double xk = (7*M_PI)/2;
double a;

double fXY(double x, double y) {
    return 5*sin(x)*( (5*cos(x)) + y );
}

double finalFX(double x) {
    return exp(-5 * cos(x)) - (5*cos(x)) + 1;
}

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    double h = atof(argv[1]);

    int mode = atoi(argv[2]);
    bool euler = false;

    if (mode == 1) {
        euler = true;
        cout << "Metoda: Euler'a\n";
    } else if (mode == 2) {
        euler = false;
        cout << "Metoda: Rungego-Kutty\n";
    } else {
        cout << "Nieprawidłowy 2 argument!\n";
        return 1;
    }
    // (xk - x0)/h  + 1 = n
    double n = ((xk-x0)/h) + 1;

    // FIXME: lepiej stała wartość
    a = finalFX(x0);

    PointDifferential * points = new PointDifferential[n];
    points[0].x = x0;
    points[0].y = a;
    points[0].dY = fXY(x0, a);

    eulerMethodDiff(points, h, n, fXY);

    sendPlotToFile(points, n, "out.dat", true);

    delete[] points;

    return 0;
}
