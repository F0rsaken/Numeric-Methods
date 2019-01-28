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
    Point * data = new Point[n];

    // punkty równoodległe
    double singleStep = (b-a)/(n-1);
    data[0].x = a;
    data[n-1].x = b;
    for (int i = 1; i < n-1; i++) {
        data[i].x = data[i-1].x + singleStep;
    }

    for (int i = 0; i < n; i++) {
        data[i].y = fX(data[i].x);
    }

    return data;
}

/** 
 * Argumenty:
 *  n - ilość punktow
 *  m - stopień wielomianu
*/

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    // int printMode = atoi(argv[1]); // 0 - cout, 1 - do pliku
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    if (m >= n) {
        cout << "Wartość m musi być mniejsza od wartości n!\n";
        return 0;
    }

    Point *data = getDataPoints(n);
    PolynomialFunction aproximatedFunc = polynomialRegression(data, n, m);
    
    int outN = ((data[n - 1].x - data[0].x) / step) + 1;
    // cout << "OutN: " << outN << endl;
    Point *approximation = new Point[outN];
    double xi = data[0].x;

    for (int i = 0; i < outN; xi += step, i++) {
        approximation[i].x = xi;
        approximation[i].y = aproximatedFunc.f(approximation[i].x);
    }

    cout << "Aproksymacja: \n";
    sendPlotToFile(approximation, outN, "out.dat", true);

    drawOriginalPlot();

    return 0;
}
