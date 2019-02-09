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
    Point *data = new Point[n];
    double singleStep = (b-a)/(n-1);
    double offset = M_PI_2;

    data[0].x = a;
    data[n-1].x = b;
    for (int i = 1; i < n-1; i++) {
        data[i].x = data[i-1].x + singleStep;
    }

    for (int i = 0; i < n; i++) {
        data[i].y = fX(data[i].x);
        data[i].x += offset; 
    }

    // double oldA = a, oldB = b;
    // a += offset; b += offset;

    // double singleStep = (b-a)/(n-1);
    // Point *data = new Point[n];
    // data[0].x = a;
    // data[n-1].x = b;
    // // cout << "a, b, step: " << data[0].x << "; " << data[n-1].x << "; " << singleStep << "\n\n";
    // for (int i = 1; i < n-1; i++) {
    //     data[i].x = data[i-1].x + singleStep;
    // }

    // for (int i = 0; i < n; i++) {
    //     data[i].y = fX(data[i].x);
    //     // cout << "x[" << i << "]: " << data[i].x << ",  " << "y[" << i << "]: " << data[i].y << endl;
    // }

    // a = oldA; b = oldB;

    sendPlotToFile(data, n, "points.dat", true);

    return data;
}


/** 
 * Argumenty:
 *  n - ilość punktow
 *  m - stopien aproksymacji
*/

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    // dataPoints musi byc przesuniete o 0.5pi
    if (n%2 == 0) { // parzysta ilosc punktow
        if (m > ((n-2)/2)) {
            // m = (n-2)/2;
            cout << "Stopien wielomianu jest za wysoki, dopuszczalny stopien: (n-2)/2 = " << (n-2)/2 << endl;
            // return 0;
        }
    } else { // nieparzysta ilosc punktow
        if ( m > ( (n-1)/2 ) ) {
            // m = (n-1)/2;
            cout << "Stopien wielomianu jest za wysoki, dopuszczalny stopien: (n-1)/2 = " << (n - 1) / 2 << endl;
            // return 0;
        }
    }

    drawOriginalPlot();
    Point *data = getDataPoints(n);

    int outN = ((b - a) / step) + 1;
    // cout << "OutN: " << outN << endl;

    double rangeM = 2.5;
    double ** abFactors = countFactorsForTrygonometricApproximation(data, n, m, rangeM);
    double xi = a;

    Point *approximation = new Point[outN];
    for (int i = 0; i < outN; xi += step, i++) {
        approximation[i].x = xi;
        approximation[i].y = trygonometricApproximation(abFactors, n, rangeM, approximation[i].x + M_PI_2);
    }
    
    cout << "Aproksymacja: \n";
    sendPlotToFile(approximation, outN, "out.dat", true);

    delete[] data;

    return 0;
}