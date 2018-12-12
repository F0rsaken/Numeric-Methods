#ifndef EQUATIONS_LIB_H
#define EQUATIONS_LIB_H

class FunctionOf3 {
public:
    double (*f)(double, double, double);
    // auto f;

    FunctionOf3() {}
    FunctionOf3( double f(double, double, double) );
    ~FunctionOf3() {};
    void setF(double f(double, double, double));
};

class FunctionOf1 {
public:
    double (*f)(double);
    // auto f;

    FunctionOf1() {}
    FunctionOf1( double f(double) );
    ~FunctionOf1() {};
    void setF(double f(double));
};

double newtonMethod(double x0, double f(double), double fPrim(double), int stopCriterium, double ro);
double eulerMethod(double a, double b, double f(double), int stopCriterium, double ro);
double * newtonEquationsSystem (FunctionOf3 fVector[], FunctionOf1 *jacobianMatrix[], double initVector[], int stopCriterium, double ro);

#endif
