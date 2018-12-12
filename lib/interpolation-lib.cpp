#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "interpolation-lib.h"
#include "matrix-lib.h"
using namespace std;

Polynomial::Polynomial(double a, double b, double c, double d) {
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
}

double Polynomial::f(double x) {
    double retVal = ( this->a * pow(x, 3) ) + (this->b * pow(x, 2)) + (this->c * x) + this->d;
    return retVal;
}

Polynomial2::Polynomial2(double a, double b, double c) {
    this->a = a;
    this->b = b;
    this->c = c;
}

double Polynomial2::f(double x) {
    double retVal = (this->a * pow(x, 2)) + (this->b * x) + this->c;
    return retVal;
}

PolynomialFunction::PolynomialFunction() {}

PolynomialFunction::PolynomialFunction(int size) {
    this->factors = new double[size];
    this->size = size;
}

PolynomialFunction::PolynomialFunction(int size, double *points) {
    this->size = size;
    this->factors = points;
}

PolynomialFunction::~PolynomialFunction() {
    delete[] this->factors;
}

double PolynomialFunction::f(double x) {
    if (this->size <= 0) {
        return 0;
    }

    double retVal = this->factors[0];
    for (int i = 1; i < this->size; i++) {
        retVal += (this->factors[i] * pow(x, i));
    }

    return retVal;
}

void sendPlotToFile(Point data[], int n, string fileName, bool informUser) {
    // cout << "Interpolacja: \n";
    ofstream outputFile;
    outputFile.open(fileName, ios::trunc);
    for(int i = 0; i < n; i++) {
        outputFile << data[i].x << " " << data[i].y << endl;
    }
    outputFile.close();
    if (informUser) cout << "Skończono pisać do pliku " << fileName << endl;
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

void fillNewtonDiffTable(int n, Point data[], double *diffTable[]) {
    for (int i = 0; i < n; i++) { diffTable[i][0] = data[i].y; }

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n-1; j++) {
            diffTable[j][i] = diffTable[j+1][i-1] - diffTable[j][i-1];
        }
    }
}


Polynomial * cubicSplines(Point *data, int n, int boundaryType) {
    // hi = xi - xi-1
    double *h = new double[n-1];
    for (int i = 0; i < n-1; i++) { h[i] = data[i+1].x - data[i].x; }

    double *lower = new double[n-1]; // mi
    double *upper = new double[n-1]; // lambda
    double tmp;
    for (int i = 0; i < n-2; i++) {
        tmp = h[i] + h[i+1];
        lower[i] = h[i]/tmp;
        upper[i+1] = h[i+1]/tmp;
    }
    double *v = new double[n]; // f[xi-1, xi, xi+1]
    for (int i = 1; i < n-1; i++) {
        v[i] = 6*(1/(h[i]+h[i-1]))*( ( (data[i+1].y - data[i].y)/h[i] ) - ( (data[i].y - data[i-1].y)/h[i-1] ) );
    }

    double *midVector = new double[n];
    for (int i = 0; i < n; i++) { midVector[i] = 2; }

    if (boundaryType == 1) { // natural boundary condition
        midVector[0] = 1; midVector[n-1] = 1;
        lower[n-2] = 0;
        upper[0] = 0;
        v[0] = 0; v[n-1] = 0;
    } else { // clamped boundary condition
        upper[0] = 1;
        lower[n-2] = 1;
        v[0] = 6*( (data[1].y - data[0].y)/(h[0]*( data[1].x - data[0].x ) ) );
        v[n-1] = 6*( (data[n-1].y - data[n-2].y)/(h[n-2]*( data[n-1].x - data[n-2].x ) ) );
    }

    double *M = thomasAlgorithm(lower, midVector, upper, v, n);
    Polynomial *retSplines = new Polynomial[n-1];

    for (int i = 0; i < n-1; i++) {
        retSplines[i].a = (M[i+1]-M[i])/(6 * h[i]);
        retSplines[i].b = ( ( (data[i+1].x * M[i]) - (data[i].x * M[i+1]) )/(2 * h[i]) );
        retSplines[i].c = (( data[i+1].y - data[i].y )/h[i]) - ( ( (h[i] * M[i+1]) - (h[i] * M[i]) )/6 );
        retSplines[i].d = (( (data[i+1].x * data[i].y) - (data[i].x * data[i+1].y) )/h[i]) - (( (h[i] * data[i+1].x * M[i]) - (h[i] * data[i].x * M[i+1]) )/6);
    }

    delete[] h;
    delete[] lower;
    delete[] upper;
    delete[] midVector;
    delete[] v;
    delete[] M;

    return retSplines;
}

Polynomial2 * quadraticSplines(Point *data, int n) {
    double *h = new double[n-1];
    for (int i = 0; i < n-1; i++) { h[i] = data[i+1].x - data[i].x; }

    double *D = new double[n];
    D[0] = 0;
    for (int i = 1; i < n; i++) {
        D[i] = ( ( (2*data[i].y) - (2*data[i-1].y) )/h[i-1] ) - D[i-1];
    }

    Polynomial2 *retSplines = new Polynomial2[n-1];

    for (int i = 0; i < n-1; i++) {
        retSplines[i].a = ( (D[i+1] - D[i]) / (2*h[i]) );
        retSplines[i].b = ( ( (data[i+1].x * D[i]) - (data[i].x * D[i-1]) ) / h[i] );
        retSplines[i].c = ( data[i+1].y - (D[i+1] * h[i] * 2) );
    }

    delete[] h;
    delete[] D;

    return retSplines;
}

PolynomialFunction polynomialRegression(Point *data, int n) {
    int xSize = (2*n) - 1;
    double *xCountedValues = new double[xSize];
    xCountedValues[0] = 1;
    // robimy tablice, aby przyspieszyć liczenie
    for (int i = 1; i < xSize; i++) {
        xCountedValues[i] = 0;
        for (int j = 0; j < n; j++) {
            xCountedValues[i] =pow(data[j].x, i);
        }
    }

    double **aMatrix = new double*[n];
    for (int i = 0; i < n; i++) { aMatrix[i] = new double[n]; }

    for (int j = 0; j < n; j++) {
        for (int i = j; i < j+n; i++) {
            aMatrix[j][i] = xCountedValues[i];
        }
    }
    aMatrix[0][0] = n;
    double *bVector = new double[n];
    // bVector[0] = 0;
    for (int i = 0; i < n; i++) {
        bVector[i] = 0;
        for (int j = 0; j < n; i++) {
            bVector[i] += (data[j].y * pow(data[j].x, (double)j));
        }
    }

    double *initVector = new double[n];
    for (int i = 0; i < n; i++) { initVector[i] = 0; }
    double *xVector = SORAlgorithm(n, (const double**)aMatrix, initVector, bVector, 0.01, 1, 1);

    PolynomialFunction aproximation = PolynomialFunction(n, xVector);

    for (int i = 0; i < n; i++) {
        delete[] aMatrix[i];
    }
    delete[] aMatrix;
    delete[] xCountedValues;
    delete[] bVector;

    return aproximation;
}

double ** countFactorsForTrygonometricApproximation(Point *dataPoints, int n) {
    double **factorsTable = new double*[2];
    factorsTable[0] = new double[n];
    factorsTable[1] = new double[n];

    for (int i = 0; i < n; i++) {
        // coś * 2 i potem przez n
        factorsTable[0][i] = 0; factorsTable[1][i] = 0;
        for (int j = 0; j <= i; j++) {
            // punkt a
            factorsTable[0][i] += ( dataPoints[j].y * cos( i*dataPoints[j].x ) );
            // punkt b
            factorsTable[1][i] += ( dataPoints[j].y * sin( i*dataPoints[j].x ) );
        }
        factorsTable[0][i] *= 2; factorsTable[1][i] *= 2;

        factorsTable[0][i] /= n; factorsTable[1][i] /= n;
    }

    return factorsTable;
}

double trygonometricApproximation(double *factorsTable[], int n, Point *dataPoints, double x) {
    double *aFactors = factorsTable[0];
    double *bFactors = factorsTable[1];
    double outVal = aFactors[0]/2;

    double aSum = 0;
    double bSum = 0;
    for (int i = 1; i < n; i++) {
        aSum += (aFactors[i] * cos(i * x));
        bSum += (bFactors[i] * sin(i * x));
    }
    outVal += aSum + bSum;

    return outVal;
}