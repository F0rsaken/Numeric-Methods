#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <math.h>
// #include "../../lib/matrix-lib.h"
#include "../../lib/interpolation-lib.h"
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
    double tmp1 = 5*sin(x);
    double tmp2 = 5*cos(x);
    return (tmp1 * (tmp2 + y));
}

double finalFX(double x) {
    double tmp1 = exp((-5) * cos(x));
    double tmp2 = (5*cos(x)) + 1;
    return (tmp1 - tmp2);
}

void drawOriginalPlot() {
    double delta = 0.05;
    int n = (xk - x0)/delta;
    Point *originalPoints = new Point[n];
    originalPoints[0].x = x0;
    originalPoints[0].y = finalFX(x0);
    for (int i = 1; i < n; i++) {
        originalPoints[i].x = originalPoints[i - 1].x + delta;
        originalPoints[i].y = finalFX(originalPoints[i].x);
    }

    sendPlotToFile(originalPoints, n, "original_plot.dat", true);

    delete[] originalPoints;
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
    int n = ((xk-x0)/h) + 1;

    // FIXME: lepiej stała wartość
    a = finalFX(x0);

    PointDifferential * points = new PointDifferential[n];
    points[0].x = x0;
    points[0].y = a;
    points[0].dY = fXY(x0, a);

    if (euler) {
        eulerMethodDiff(points, h, n, fXY);
    } else {
        rungeKuttyMethodDiff(points, h, n, fXY);
    }

    // TODO: ograniczyc ilosc punktow dla wyjścia
    sendPlotToFileDiff(points, n, "out.dat", true);

    drawOriginalPlot();

    delete[] points;

    return 0;
}
