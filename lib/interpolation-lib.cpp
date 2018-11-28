#include <cmath>
#include "interpolation-lib.h"
#include "matrix-lib.h"

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