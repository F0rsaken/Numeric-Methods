#ifndef MATRIX_LIB_H
#define MATRIX_LIB_H

// lab 1
void printMatrix(int m, int n, double** matrix);
void printVector(int n, const double *vector, bool vertical = true);
double vectorEuclideanNorm(const double * vector, int length);
double * generateXVector(int n, const double dataTab[], int dataN);
double * multiplyMatrixByVector(int m, int n, const double *matrix[], double vector[]);
double * multiplyTridiagonalMatrixByVector(int n, double *matrix[], double vector[]);
double ** multiplyMatrixByMatrix(int n, double *matrixA[], double *matrixB[]);
double * gaussElimination(int rows, int columns, double *matrix[], double vector[]);
double * gaussEliminationOld(int rows, int columns, double *matrix[], double vector[]);
double * thomasAlgorithm(const double *lower, const double *mid, const double *upper, const double *freeVector, int N);
double * thomasAlgorithmOld(int columns, double *matrix[], double vector[]);

//
float * generateXVector_f(int n, const float dataTab[], int dataN);
float * multiplyMatrixByVector_f(int m, int n, float *matrix[], float vector[]);
float * gaussElimination_f(int rows, int columns, float *matrix[], float vector[]);
float vectorEuclideanNorm_f(float * vector, int length);
//

// lab 2
// a)
double * jacobiAlgorithm(int n, const double *matrix[], double *initVectorX, double *vectorB, double ro, int stopType);
double * substractVectors(const double *vector1, const double *vector2, int n);
// b)
double powerIteration(double *matrix[], double *vector, int n, double stopVal);
double rayleighQuotient(double *matrix[], double *vector, int n);
// c)
double * SORAlgorithm(int n, const double *matrix[], double *initVectorX, double *vectorB, double ro, double relax, int stopType);

#endif