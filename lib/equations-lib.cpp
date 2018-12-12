#include <cmath>
#include <iostream>
#include "equations-lib.h"
#include "matrix-lib.h"
using namespace std;

FunctionOf3::FunctionOf3(double f(double, double, double)) {
    this->f = f;
}

void FunctionOf3::setF(double f(double, double, double)) {
    this->f = f;
}

FunctionOf1::FunctionOf1(double f(double)) {
    this->f = f;
}

void FunctionOf1::setF(double f(double)) {
    this->f = f;
}

double newtonMethod (double x0, double f(double), double fPrim(double), int stopCriterium, double ro) {
    bool stop = false;
    double y, yPrim, xk;
    int i;

    for (i = 0; !stop; i++) {
        y = f(x0);
        yPrim = fPrim(x0);

        xk = x0 - (y/yPrim);
        if (stopCriterium == 1) {
            if ( fabs(xk - x0) < ro ) {
                stop = true;
            }
        } else {
            if ( fabs(f(x0)) < ro ) {
                stop = true;
            }
        }

        if (!stop) {
            x0 = xk;
            xk = 0;
        }
    }

    cout << "Przeprowadzono " << "\033[35m" << i << "\033[0m" << " iteracji\n";

    return xk;
}

double eulerMethod(double a, double b, double f(double), int stopCriterium, double ro) {
    double x0 = a, x1 = b, xn, counter1, counter2, denominator;
    bool stop = false;
    int i;
    for (i = 0; !stop; i++) {
        counter1 = f(x1) * x0;
        counter2 = f(x0) * x1;
        denominator = f(x1) - f(x0);
        xn = (counter1/denominator) - (counter2/denominator);

        if (stopCriterium == 1) {
            if ( fabs(xn - x0) < ro ) {
                stop = true;
            }
        } else {
            if ( fabs(f(x0)) < ro ) {
                stop = true;
            }
        }

        if (!stop) {
            x0 = x1;
            x1 = xn;
            xn = 0;
        }
    }

    cout << "Przeprowadzono " << "\033[35m" << i << "\033[0m" << " iteracji\n";

    return xn;
}

double *newtonEquationsSystem(FunctionOf3 fVector[], FunctionOf1 *jacobianMatrix[], double initVector[], int stopCriterium, double ro) {
    double *fVectorValues = new double[3];
    double **jacobianMatrixValues = new double*[3];
    double *initVectorCopy = new double[3];
    for (int i = 0; i < 3; i++) { jacobianMatrixValues[i] = new double[3]; initVectorCopy[i] = initVector[i]; }

    auto countFVector = [fVector](double *fVal, double *vec) {
        for (int i = 0; i < 3; i++) {
            fVal[i] = 0;
            fVal[i] = fVector[i].f(vec[0], vec[1], vec[2]);
        }
    };
    auto countJacobian = [jacobianMatrix](double **matrix, double *vec) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                matrix[i][j] = 0;
                matrix[i][j] = jacobianMatrix[i][j].f(vec[j]);
                // cout << "x: " << vec[j] << ", f(x): " << matrix[i][j] << endl;
            }
        }
    };

    bool stop = false;
    double * Xn, *tmpVector, SORInitVector[3] = {0, 0, 0};
    for (int i = 0; !stop; i++) {
        cout << "Iter number: " << i << endl;
        countFVector(fVectorValues, initVectorCopy);
        for (int j = 0; j < 3; j++) { fVectorValues[j] *= -1; }
        // printVector(3, fVectorValues);
        countJacobian(jacobianMatrixValues, initVectorCopy);
        // printMatrix(3, 3, jacobianMatrixValues);
        Xn = gaussElimination(3, 3, jacobianMatrixValues, fVectorValues);
        // Xn = jacobiAlgorithm(3, (const double**)jacobianMatrixValues, SORInitVector, fVectorValues, 0.1, 1);
        for (int j = 0; j < 3; j++) { Xn[j] += initVectorCopy[j]; }
        // printVector(3, Xn);

        if (stopCriterium == 1) {
            if ( vectorEuclideanNorm( tmpVector = substractVectors(Xn, initVectorCopy, 3), 3 ) < ro ) {
                stop = true;
            } else {
                cout << "Still counting..." << endl;
            }
            cout << vectorEuclideanNorm( tmpVector, 3 ) << endl;
            delete[] tmpVector;
        } else {
            if ( vectorEuclideanNorm(fVectorValues, 3) < ro ) {
                stop = true;
            } else {
                cout << "Still counting..." << endl;
            }
        }

        if (!stop) {
            tmpVector = initVectorCopy;
            initVectorCopy = Xn;
            // Xn = NULL;
            // cout << "Init vector:\n";
            // printVector(3, initVectorCopy);
            delete[] tmpVector;
            // cout << "Init vector:\n";
            // printVector(3, initVectorCopy);
            // if (Xn == NULL) cout << "To null, czyli gut \n\n";
        }
        // Xn = SORAlgorithm(3, jacobianMatrixValues, tmpVector, fVectorValues, 0.0001, 1, 1);
    }

    delete[] initVectorCopy;
    delete[] fVectorValues;
    for (int i = 0; i < 3; i++) { delete[] jacobianMatrixValues[i]; }
    delete[] jacobianMatrixValues;

    return Xn;
}
// reset             0  (everything back to normal)
// bold/bright       1  (often a brighter shade of the same colour)
// underline         4
// inverse           7  (swap foreground and background colours)
// bold/bright off  21
// underline off    24
// inverse off      27

//          foreground background
// black        30         40
// red          31         41
// green        32         42
// yellow       33         43
// blue         34         44
// magenta      35         45
// cyan         36         46
// white        37         47