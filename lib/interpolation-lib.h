#ifndef INTERPOL_LIB_H
#define INTERPOL_LIB_H
#include <string>
using namespace std;

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


class Polynomial {
public:
    double a, b, c, d, x;

    Polynomial() {};
    Polynomial(double a, double b, double c, double d);
    ~Polynomial() {};

    double f(double x);
};


class Polynomial2 {
public:
    double a, b, c, x;

    Polynomial2(){};
    Polynomial2(double a, double b, double c);
    ~Polynomial2(){};

    double f(double x);
};

class PolynomialFunction {
public:
    double *factors;
    int size;

    PolynomialFunction();
    PolynomialFunction(int size);
    PolynomialFunction(int size, double *points);
    ~PolynomialFunction();

    double f(double x);
};

void sendPlotToFile(Point data[], int n, string fileName, bool informUser = false);

// functions
double useLagrange(Point data[], int n, double xi);
double useNewton(Point data[], int n, double x, double *diffTable[]);
void fillNewtonDiffTable(int n, Point data[], double *diffTable[]);
Polynomial *cubicSplines(Point *data, int n, int boundaryType);
Polynomial2 * quadraticSplines(Point *data, int n);
PolynomialFunction polynomialRegression(Point *data, int n);

// aproksymacje
double ** countFactorsForTrygonometricApproximation(Point * dataPoints, int n);
double trygonometricApproximation(double *factorsTable[], int n, Point * dataPoints, double x);

#endif