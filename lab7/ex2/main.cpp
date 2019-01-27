#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include <math.h>
#include "../../lib/matrix-lib.h"
#include "../../lib/interpolation-lib.h"
#include "../../lib/diff-equations-lib.h"
using namespace std;


/**
 * FIXME: inna funkcja
 * Zadana funkcja:
 * y'' + 4y = 4x
 * y(0) = 0
 * y((2pi+1)/2) = ?
 * 
 * Zadane parametry:
 * x0 = 0
 * xk = (2PI + 1)/2
 * 
 * 
 * rozwiazanie:
 *  y(x) = -sin(2x) + x;
 */

/**
 * Argumenty:
 *  - h
 **/

double xS = 0;
double xE = ((2*M_PI)+1)/2;
double yS = 0;
double yE;
double step = 0.02;

// TODO:
double finalFX(double x) {
    double tmp1 = -1 * sin(2*x);
    return (tmp1 + x);
}

void drawOriginalPlot() {
    int n = ((xE - xS) / step) + 1;
    Point *dataPoints = new Point[n];
    double x = xS;
    for (int i = 0; i < n-1; i++, x += step) {
        dataPoints[i].x = x;
        dataPoints[i].y = finalFX(x);
    }
    dataPoints[n-1].x = xE;
    dataPoints[n-1].y = finalFX(xE);
    sendPlotToFile(dataPoints, n, "original_plot.dat", true);
    delete[] dataPoints;
}

void getBestCountedPoints(int n, Point points[]) {
    int outN = ((xE - xS) / step) + 1;
    if (outN > n) {
        sendPlotToFile(points, n, "out.dat", true);
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
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    yE = finalFX(xE);
    double h = atof(argv[1]);
    cout << "1/h^2: " << 1/pow(h,2) << endl;

    // (xk - x0)/h  + 1 = n
    int n = ((xE-xS)/h) + 1;

    Point startPoint = Point(xS, yS);
    Point endPoint = Point(xE, yE);

    Point *points = MRSAlgorithm(h, n, startPoint, endPoint);

    getBestCountedPoints(n, points);

    drawOriginalPlot();

    delete[] points;

    return 0;
}
