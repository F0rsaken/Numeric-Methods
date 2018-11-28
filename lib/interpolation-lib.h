#ifndef INTERPOL_LIB_H
#define INTERPOL_LIB_H

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

// functions
Polynomial * cubicSplines(Point *data, int n, int boundaryType);
Polynomial2 * quadraticSplines(Point *data, int n);

#endif