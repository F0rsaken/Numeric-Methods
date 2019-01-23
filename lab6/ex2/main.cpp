#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <tgmath.h>
#include "../../lib/matrix-lib.h"
#include "../../lib/interpolation-lib.h"
using namespace std;

double step = 0.1;

double a = -3 * M_PI;
double b = 2 * M_PI;

// y = e^(-3( sin(x) ))
double fX(double x) {
    return exp( -3 * sin(x) );
}

void drawOriginalPlot() {
    int n = ((b - a) / step) + 1;
    Point *dataPoints = new Point[n];
    double x = a;
    for (int i = 0; i < n-1; i++, x += step) {
        dataPoints[i].x = x;
        dataPoints[i].y = fX(x);
    }
    dataPoints[n-1].x = b;
    dataPoints[n-1].y = fX(b);
    sendPlotToFile(dataPoints, n, "original_plot.dat", true);
}


Point* getDataPoints(int n) {
    // przesuniecie funkcji
    double oldA = a, oldB = b;
    a = 0; b = 2 * M_PI;

    double singleStep = (b-a)/(n-1);
    Point *data = new Point[n];
    data[0].x = a;
    data[n-1].x = b;
    // cout << "a, b, step: " << data[0].x << "; " << data[n-1].x << "; " << singleStep << "\n\n";
    for (int i = 1; i < n-1; i++) {
        data[i].x = data[i-1].x + singleStep;
    }

    for (int i = 0; i < n; i++) {
        data[i].y = fX(data[i].x);
        // cout << "x[" << i << "]: " << data[i].x << ",  " << "y[" << i << "]: " << data[i].y << endl;
    }

    a = oldA; b = oldB;

    return data;
}

/** 
 * Argumenty:
 *  n - ilość punktow
*/

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    int n = atoi(argv[1]);

    drawOriginalPlot();

    Point *data = getDataPoints(n);
    int outN = ((b - a) / step) + 1;
    // cout << "OutN: " << outN << endl;

    Point *approximation = new Point[outN];
    double ** abFactors = countFactorsForTrygonometricApproximation(data, n);
    double xi = a;

    for (int i = 0; i < outN; xi += step, i++) {
        approximation[i].x = xi;
        approximation[i].y = trygonometricApproximation(abFactors, n, data, approximation[i].x);
    }
    
    cout << "Aproksymacja: \n";
    sendPlotToFile(approximation, outN, "out.dat", true);

    // ofstream outputFile;
    // outputFile.open("out.dat", ios::trunc);
    // for(int i = 0; i < outN; i++) {
    //     outputFile << approximation[i].x << " " << approximation[i].y << endl;
    //     // cout << "x: " << interpolation[i].x << ", y: " << interpolation[i].y << endl;
    // }
    // outputFile.close();
    // cout << "Skończono pisać do pliku\n";

    delete[] data;

    return 0;
}
