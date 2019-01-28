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
    double a, b, c, d;

    Polynomial() {};
    Polynomial(double a, double b, double c, double d);
    ~Polynomial() {};

    double f(double x);
};

class CubicSpline {
  public:
    double h, x1, x2, y1, y2, M1, M2;

    CubicSpline() {};
    CubicSpline(double h, double x1, double x2, double y1, double y2, double M1, double M2) {
        this->h = h;
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
        this->M1 = M1;
        this->M2 = M2;
    }
    ~CubicSpline() {};

    double fx(double x);
};

class QuadraticSpline {
  public:
    double h, x1, x2, y, D1, D2;

    QuadraticSpline() {};
    QuadraticSpline(double h, double x1, double x2, double y, double D1, double D2) {
        this->h = h;
        this->x1 = x1;
        this->x2 = x2;
        this->y = y;
        this->D1 = D1;
        this->D2 = D2;
    }
    ~QuadraticSpline() {};

    double fx(double x);

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
double useHermit(Point data[], int n, double x, double *diffTable[]);
void fillNewtonDiffTable(int n, Point data[], double *diffTable[]);
void fillHermitDiffTable(int n, Point data[], double *diffTable[], double derivatives[]);
CubicSpline *cubicSplines(Point *data, int n, int boundaryType);
QuadraticSpline * quadraticSplines(Point *data, int n, double D0);

// aproksymacje
PolynomialFunction polynomialRegression(Point *data, int n, int m);
double ** countFactorsForTrygonometricApproximation(Point * dataPoints, int nPoints, int degree);
double trygonometricApproximation(double *factorsTable[], int degree, double x);

#endif