#include <cstdlib>
#include <iostream>
#include <time.h>
#include <array>
#define _USE_MATH_DEFINES
#include <math.h>
#include <tgmath.h>
#include <fstream>
#include <string>
#include "../../lib/interpolation-lib.h"
#include "../../lib/matrix-lib.h"
using namespace std;

double step = 0.2;

double a = -3 * M_PI;
double b = 2 * M_PI;

// y = e^(-3( sin(x) ))
double fX(double x) {
    return exp( -3 * sin(x) );
}
double fPrimX(double x) {
    return -3 * exp(-3 * sin(x)) * cos(x);
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

Point* getDataPoints(Point data[], int n, bool czebyszew) {
    data = new Point[n];
    if (czebyszew) {
        // zera wielomianu Czybyszewa
        for (int i = 0; i < n; i++) {
            data[n-i-1].x = ( ((b - a)/2) * cos( ((2*i + 1) * M_PI)/(2 * n) ) ) + ((b + a)/2);
        }

        data[0].x = a;
        data[n-1].x = b;
    } else {
        // punkty równoodległe
        double singleStep = (b-a)/(n-1);
        data[0].x = a;
        data[n-1].x = b;
        for (int i = 1; i < n-1; i++) {
            data[i].x = data[i-1].x + singleStep;
        }
    }
    for (int i = 0; i < n; i++) {
        data[i].y = fX(data[i].x);
    }

    return data;
}

double* getDerivatives(Point data[], int n) {
    double * retTab = new double[n];

    for (int i = 0; i < n; i++) {
        retTab[i] = fPrimX(data[i].x);
    }

    return retTab;
}

/** 
 * Argumenty:
 *  0 | 1 - rozłożenie punktów: równoodległe | zera Czybyszewa
 *  n - ilość punktow
*/

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    // parsowanie argumentów
    int pointsMode = atoi(argv[1]); // 0 - rownoodlegle, 1 - zera Czebyszewa
    int n = atoi(argv[2]);

    bool lagrange = false;
    bool czebyszew = false;

    if (pointsMode == 0) {
        czebyszew = false;
    } else if (pointsMode == 1) {
        czebyszew = true;
    } else {
        cout << "Zła wartość pierwszego argumentu!\n";
        return 1;
    }

    // generowanie punktów
    Point *data = getDataPoints(data, n, czebyszew);

    cout << "Data points:\n";
    for (int i = 0; i < n; i++) {
        cout << "x: " << data[i].x << ", y: " << data[i].y << endl;
    }

    double *derivatives = getDerivatives(data, n);

    // tabela x do przedstawienia wykresu do interpolacji
    double xi = data[0].x;
    int outN = fabs((data[n-1].x - data[0].x)/step) + 1;
    cout << "OutN: " << outN << endl;
    Point *interpolation = new Point[outN];

    double **forwardDiffTable;

    // stwórz tablice
    forwardDiffTable = new double*[2*n];
    for (int i = 0; i < 2*n; i++) { forwardDiffTable[i] = new double[(2*n)+1]; }
    // // wypełnij tablice
    fillHermitDiffTable(n, data, forwardDiffTable, derivatives);

    for (int i = 0; i < outN; xi += step, i++) {
        interpolation[i].x = xi;
        interpolation[i].y = useHermit(data, n, xi, forwardDiffTable);
    }
    interpolation[outN-1].x = data[n-1].x;
    interpolation[outN-1].y = useHermit(data, n, b, forwardDiffTable);


    drawOriginalPlot();

    cout << "Interpolacja: \n";
    sendPlotToFile(interpolation, outN, "out.dat", true);

    return 0;
}