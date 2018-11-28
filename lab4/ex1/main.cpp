#include <cstdlib>
#include <iostream>
#include <time.h>
#include <array>
#define _USE_MATH_DEFINES
#include <math.h>
#include <tgmath.h>
#include <fstream>
// #include "../../lib/matrix-lib.h"
using namespace std;

double step = 0.1;

double a = -3 * M_PI;
double b = 2 * M_PI;

class Point {
public:
    double x;
    double y;

    Point() {}
    Point(double x, double y) {
        this->x = x;
        this->y = y;
    }
};

// y = e^(-3( sin(x) ))
double fX(double x) {
    return exp( -3 * sin(x) );
}

Point* getDataPoints(Point data[], int n) {
    double singleStep = (b-a)/(n-1);
    data = new Point[n];
    data[0].x = a;
    data[n-1].x = b;
    // cout << "a, b, step: " << data[0].x << "; " << data[n-1].x << "; " << singleStep << "\n\n";
    for (int i = 1; i < n-1; i++) {
        data[i].x = data[i-1].x + singleStep;
    }

    for (int i = 0; i < n; i++) {
        data[i].y = fX(data[i].x);
        cout << "x[" << i << "]: " << data[i].x << ",  " << "y[" << i << "]: " << data[i].y << endl;;
    }

    return data;
}

double useLagrange(Point data[], int n, double xi) {
    double retVal = 0, upper, lower, li;
    for (int i = 0; i < n; i++) {

        upper = 1, lower = 1;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                upper = upper * (xi - data[j].x);
                lower = lower * (data[i].x - data[j].x);
            }
        }
        li = upper/lower;
        retVal += (data[i].y * li);
    }
    return retVal;
}

void fillNewtonDiffTable(int n, Point data[], double *diffTable[]) {
    for (int i = 0; i < n; i++) { diffTable[i][0] = data[i].y; }

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n-1; j++) {
            diffTable[j][i] = diffTable[j+1][i-1] - diffTable[j][i-1];
        }
    }
}

double useNewton(Point data[], int n, double x, double *diffTable[]) {
    double retVal = diffTable[0][0];
    double ni;

    for (int i = 1; i < n; i++) {
        ni = 1;
        for (int j = 0; j < i; j++) {
            ni = ni * (x-data[j].x);
        }
        retVal += diffTable[0][i] * ni;
    }

    return retVal;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    // int printMode = atoi(argv[1]); // 0 - cout, 1 - do pliku
    int calcMode = atoi(argv[1]); // 0 - lagrange, 1 - Newton
    int n = atoi(argv[2]);
    bool lagrange = false;
    
    if (calcMode == 0) {
        lagrange = true;
    } else if (calcMode == 1) {
        lagrange = false;
    } else {
        cout << "Zła wartość pierwszego argumentu!\n";
        return 1;
    }

    Point *data = getDataPoints(data, n);

    cout << "Data points:\n";
    for (int i = 0; i < n; i++) {
        cout << "x: " << data[i].x << ", y: " << data[i].y << endl;
    }

    double xi = data[0].x;
    int outN = ((data[n-1].x - data[0].x)/step) + 1;
    cout << "OutN: " << outN << endl;
    Point *interpolation = new Point[outN];
    double **forwardDiffTable;

    if (lagrange) {
        for (int i = 0; i < outN; xi += step, i++) {
            interpolation[i].x = xi;
            interpolation[i].y = useLagrange(data, n, xi);
        }
    } else {
        forwardDiffTable = new double*[n];
        for (int i = 0; i < n; i++) { forwardDiffTable[i] = new double[n]; }
        fillNewtonDiffTable(n, data, forwardDiffTable);

        for (int i = 0; i < outN; xi += step, i++) {
            interpolation[i].x = xi;
            interpolation[i].y = useNewton(data, n, xi, forwardDiffTable);
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

    return 0;
}