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
 * FIXME: inna funkcja
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
 * a = 2
 * 
 * rozwiazanie:
 *  y(x) = e^(-5 * cos(x)) - 5*cos(x) + 1
 */

/**
 * Argumenty:
 *  - h
 **/

double x0 = M_PI/2;
double xk = (7*M_PI)/2;
double a;
double step = 0.05;

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
    int n = ((xk - x0) / step) + 1;
    Point *dataPoints = new Point[n];
    double x = x0;
    for (int i = 0; i < n-1; i++, x += step) {
        dataPoints[i].x = x;
        dataPoints[i].y = finalFX(x);
    }
    dataPoints[n-1].x = xk;
    dataPoints[n-1].y = finalFX(xk);
    sendPlotToFile(dataPoints, n, "original_plot.dat", true);
    delete[] dataPoints;
}

void getBestCountedPoints(int n, PointDifferential points[]) {
    int outN = ((xk - x0) / step) + 1;
    if (outN > n) {
        sendPlotToFileDiff(points, n, "out.dat", true);
        return;
    }

    Point *outPoints = new Point[outN];

    outPoints[0].x = points[0].x;
    outPoints[0].y = points[0].y;

    double x = outPoints[0].x;
    for (int i = 1, j = 1; i < outN-1; i++) {
        x += step;
        while(x > points[j].x ) { j++; }

        outPoints[i].x = points[j].x;
        outPoints[i].y = points[j].y;
    }

    outPoints[outN-1].x = points[n-1].x;
    outPoints[outN-1].y = points[n-1].y;
    sendPlotToFile(outPoints, outN, "out.dat", true);
}

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    double h = atof(argv[1]);

    // (xk - x0)/h  + 1 = n
    int n = ((xk-x0)/h) + 1;

    // FIXME: lepiej stała wartość
    a = 2;

    PointDifferential * points = new PointDifferential[n];
    points[0].x = x0;
    points[0].y = a;
    points[0].dY = fXY(x0, a);

    // FIXME: rozwiazanie

    getBestCountedPoints(n, points);

    drawOriginalPlot();

    delete[] points;

    return 0;
}
