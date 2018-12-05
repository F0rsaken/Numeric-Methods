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

double step = 0.01;

double a = -3 * M_PI;
double b = 2 * M_PI;

// y = e^(-3( sin(x) ))
double fX(double x) {
    return exp( -3 * sin(x) );
}

Point* getDataPoints(int n) {
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

    return data;
}

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    // int printMode = atoi(argv[1]); // 0 - cout, 1 - do pliku
    int boundaryType = atoi(argv[1]); // 0 - clamped, 1 - natural
    int n = atoi(argv[2]);
    
    if (boundaryType != 0 && boundaryType != 1) {
        cout << "Zła wartość pierwszego argumentu!\n";
        return 1;
    }

    Point *data = getDataPoints(n);
    Polynomial *splines = cubicSplines(data, n, boundaryType);
    
    int outN = ((data[n - 1].x - data[0].x) / step) + 1;
    // cout << "OutN: " << outN << endl;
    Point *interpolation = new Point[outN];
    double xi = data[0].x;

    for (int i = 0; i < outN; xi += step, i++) {
            interpolation[i].x = xi;

            for (int j = 1; j < n; j++) {
                if ( xi < data[j].x ) {
                    interpolation[i].y = splines[j-1].f(xi);
                }
            }
        }
    
    cout << "Interpolacja: \n";
    ofstream outputFile;
    outputFile.open("out.dat", ios::trunc);
    for(int i = 0; i < outN; i++) {
        outputFile << interpolation[i].x << " " << interpolation[i].y << endl;
        // cout << "x: " << interpolation[i].x << ", y: " << interpolation[i].y << endl;
    }
    outputFile.close();
    cout << "Skończono pisać do pliku\n";

    delete[] data;
    delete[] splines;

    return 0;
}
