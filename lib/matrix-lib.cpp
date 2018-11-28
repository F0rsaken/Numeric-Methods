#include "matrix-lib.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

/*---------------------------------------------------------------*/
/*----------------------------LAB 1------------------------------*/

void swapRows(int start, int end, double *matrix[], double *vector) {
    double *tmp = matrix[end];
    matrix[end] = matrix[start];
    matrix[start] = tmp;
    double tmpVal = vector[end];
    vector[end] = vector[start];
    vector[start] = tmpVal;
}

void swapRows_f(int start, int end, float *matrix[], float *vector) {
    float *tmp = matrix[end];
    matrix[end] = matrix[start];
    matrix[start] = tmp;
    float tmpVal = vector[end];
    vector[end] = vector[start];
    vector[start] = tmpVal;
}

void printMatrix(int m, int n, double **matrix) {
    for (int i = 0; i < m; i++) {
        cout << "\t";
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printVector(int n, const double *vector, bool vertical) {
    if (vertical) {
        for (int i = 0; i < n; i++) {
            cout << "\t" << vector[i] << "\n";
        }
    } else {
        for (int i = 0; i < n; i++) {
            cout << "\t" << vector[i] << " ";
        }
        cout << endl;
    }
}

float * generateXVector_f(int n, const float dataTab[], int dataN) {
    float *outputMatrix = new float[n];

    int j;
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        j = rand() % dataN;
        outputMatrix[i] = dataTab[j];
    }
    return outputMatrix;
}

float * multiplyMatrixByVector_f(int rowN, int columnN, float *matrix[], float vector[]) {
    float *outVector = new float[columnN];

    for (int i = 0; i < rowN; i++) {
        outVector[i] = 0;
        for (int j = 0; j < columnN; j++) {
            outVector[i] += (matrix[i][j] * vector[j]);
        }
    }

    return outVector;
}

float * gaussElimination_f(int rows, int columns, float *matrix[], float vector[]) {
    int n = rows;
    float *solution = new float[rows];

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        float maxEl = abs(matrix[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(matrix[k][i]) > maxEl) {
                maxEl = abs(matrix[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        swapRows_f(maxRow, i, matrix, vector);

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            float c = -matrix[k][i]/matrix[i][i];
            for (int j=i; j<n; j++) {
                if (i==j) {
                    matrix[k][j] = 0;
                } else {
                    matrix[k][j] += c * matrix[i][j];
                }
            }
            vector[k] += c * vector[i];
        }
    }

    float tmp;
    for (int i = rows-1; i >= 0; i--) {
        tmp = vector[i];
        for (int j = rows-1; j > i; j--) {
            tmp -= matrix[i][j] * solution[j];
        }
        solution[i] = tmp / matrix[i][i];
    }
    return solution;
}

float vectorEuclideanNorm_f(const float * vector, int length) {
    float retVal = 0;
    for (int i = 0; i < length; i++) {
        retVal += pow(vector[i], 2);
    }
    return sqrt(retVal);
}

double vectorEuclideanNorm(const double * vector, int length) {
    double retVal = 0;
    for (int i = 0; i < length; i++) {
        retVal += pow(vector[i], 2);
    }
    return sqrt(retVal);
}

double *generateXVector( int n, const double dataTab[], int dataN ) {
    double *outputMatrix = new double[n];

    int j;
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        j = rand() % dataN;
        outputMatrix[i] = dataTab[j];
    }
    return outputMatrix;
}

double *multiplyMatrixByVector(int rowN, int columnN, const double *matrix[], double vector[]) {
    double *outVector = new double[columnN];

    for (int i = 0; i < rowN; i++) {
        outVector[i] = 0;
        for (int j = 0; j < columnN; j++) {
            outVector[i] += (matrix[i][j] * vector[j]);
        }
    }

    return outVector;
}

double *multiplyTridiagonalMatrixByVector(int n, double *matrix[], double vector[]) {
    double *outVector = new double[n];

    for (int i = 0; i<n; i++) {
        outVector[i] = matrix[1][i] * vector[i];
        if ( i < n-1 ) outVector[i] += matrix[0][i] * vector[i+1];
        if ( i > 0 ) outVector[i] += matrix[2][i] * vector[i-1];
    }

    return outVector;
}

double **multiplyMatrixByMatrix(int n, double *matrixA[], double *matrixB[]) {
    double **outMatrix = new double*[n];
    for (int i = 0; i < n; i++) { outMatrix[i] = new double[n]; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            outMatrix[i][j] = 0;
            for (int z = 0; z < n; z++) {
                outMatrix[i][j] += matrixA[i][z] * matrixB[z][j];
            }
        }
    }
    return outMatrix;
}

double *gaussEliminationOld(int rows, int columns, double *matrix[], double vector[]) {
    int pivot;
    double pivotVal, ratio;
    double *solution = new double[rows];
    // funckja nadpisuje macierz!!!!!

    for (int d = 0; d < rows-1; d++) {
        // select pivot
        pivot = d; pivotVal = matrix[pivot][d];
        for (int row = d+1; row < rows; row++) {
            if ( fabs(matrix[row][d]) > fabs(pivotVal)) {
                pivot = row;
                pivotVal = fabs(matrix[row][d]);
            }
        }

        if (pivot != d) {
            swapRows(pivot, d, matrix, vector);
        }

        for (int i = d+1; i < rows; i++) {
            ratio = matrix[i][d]/pivotVal;
            vector[i] -= vector[pivot] * ratio;

            for (int j = columns-1; j >= d; j--) {
                matrix[i][j] -= matrix[pivot][j] * ratio;
            }
            matrix[i][d] = 0;
            ratio = 0;
        }
    }

    double tmp;
    // cout << "test1\n";
    for (int i = rows-1; i >= 0; i--) {
        tmp = vector[i];
        for (int j = rows-1; j > i; j--) {
            tmp -= matrix[i][j] * solution[j];
        }
        solution[i] = tmp / matrix[i][i];
    }
    return solution;
}

double *gaussElimination(int rows, int columns, double *matrix[], double vector[]) {
    int n = rows;
    double *solution = new double[rows];

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(matrix[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(matrix[k][i]) > maxEl) {
                maxEl = abs(matrix[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        swapRows(maxRow, i, matrix, vector);

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -matrix[k][i]/matrix[i][i];
            for (int j=i; j<n; j++) {
                if (i==j) {
                    matrix[k][j] = 0;
                } else {
                    matrix[k][j] += c * matrix[i][j];
                }
            }
            vector[k] += c * vector[i];
        }
    }

    double tmp;
    for (int i = rows-1; i >= 0; i--) {
        tmp = vector[i];
        for (int j = rows-1; j > i; j--) {
            tmp -= matrix[i][j] * solution[j];
        }
        solution[i] = tmp / matrix[i][i];
    }
    return solution;
}

double * thomasAlgorithmOld(int columns, double *matrix[], double vector[]) {
    double *solution = new double[columns];
    double *upperVector = new double[columns], *freeVector = new double[columns];
    double ratio;

    upperVector[0] = matrix[0][0] / matrix[1][0];
    freeVector[0] = vector[0] / matrix[1][0];

    for (int i = 1; i < columns; i++) {
        // ratio = matrix[2][i]/matrix[1][i-1];
        ratio = 1.0 / (matrix[1][i] - (matrix[2][i] * upperVector[i-1]));
        // matrix[1][i] -= ratio * matrix[0][i-1];
        upperVector[i] = matrix[0][i] * ratio;
        // vector[i] -= ratio * vector[i-1];
        freeVector[i] = (vector[i] - (matrix[2][i] * freeVector[i-1])) * ratio;
    }

    // Backward substitution
    // solution[columns-1] = vector[columns-1]/matrix[1][columns-1];
    // for (int i = columns-2; i >= 0; i--) {
    //     solution[i] = (vector[i] - (matrix[0][i] * solution[i+1]))/matrix[1][i];
    // }
    solution[columns-1] = freeVector[columns-1];
    for (int i = columns-2; i >= 0; i--) {
        solution[i] = freeVector[i] - (upperVector[i] * freeVector[i+1]);
    }
    return solution;
}

double * thomasAlgorithm(const double *a, const double *b, const double *c, const double *d, int N) {
    double *cStar = new double[N];
    double *dStar = new double[N];

    // Pierwsza kolumna
    cStar[0] = c[0] / b[0];
    dStar[0] = d[0] / b[0];

    // Wypełnij cStar i dStar
    for (int i = 1; i < N; i++) {
        double m = 1.0 / (b[i] - a[i] * cStar[i-1]);
        cStar[i] = c[i] * m;
        dStar[i] = (d[i] - a[i] * dStar[i-1]) * m;
    }

    double *f = new double[N];
    f[N-1] = dStar[N-1];
    // cout << "f: " << endl; printVector(N, f);
    // Reverse sweep
    for (int i = N-2; i >= 0; i--) {
        f[i] = dStar[i] - (cStar[i] * f[i+1]);
    }
    // cout << "f: " << endl; printVector(N, f);
    delete[] cStar;
    delete[] dStar;

    return f;
}

/*---------------------------------------------------------------*/
/*----------------------------LAB 2------------------------------*/
double * substractVectors(const double *vector1, const double *vector2, int n) {
    double *outVector = new double[n];

    for (int i = 0; i < n; i++) {
        outVector[i] = vector1[i] - vector2[i];
    }
    return outVector;
}

bool firstStopCriterium(double ro, const double *firstVector, const double *secondVector, int n) {
    double *tmpVector = substractVectors(firstVector, secondVector, n);
    double normVal = vectorEuclideanNorm( tmpVector, n);
    delete[] tmpVector;
    return normVal < ro;
}

bool secondStopCriterium(double ro, double *vectorX, const double *vectorB, const double *matrix[], int n) {
    double *tmpVector = substractVectors( multiplyMatrixByVector(n, n, matrix, vectorX), vectorB, n);
    double normVal = vectorEuclideanNorm( tmpVector, n );
    delete[] tmpVector;
    return normVal < ro;
}

double *jacobiAlgorithm(int n, const double *matrix[], double *initVectorX, double *vectorB, double ro, int stopType) {
    bool finish = false;
    double *nextVector = new double[n];
    double sum, b, *tmpSwitch;
    int iter;
    for (iter = 1; !finish; iter++) {

        for (int i = 0; i < n; i++) {
            // sumowanie wiersz * wektor X
            sum = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += matrix[i][j] * initVectorX[j];
                }
            }

            b = 0;
            b = vectorB[i]/matrix[i][i];
            sum = sum/matrix[i][i];
            nextVector[i] = b - sum;
        }

        if (stopType == 1) {
            finish = firstStopCriterium(ro, initVectorX, nextVector, n);
        } else {
            finish = secondStopCriterium(ro, nextVector, vectorB, matrix, n);
        }

        if (!finish) {
            tmpSwitch = nextVector;
            nextVector = initVectorX;
            initVectorX = tmpSwitch;
            for (int i = 0; i < n; i++) nextVector[i] = 0;
        }

    }

    cout << "Przeprowadzono " << iter << " iteracji" << endl;

    return nextVector;
}

double powerIteration(double *matrix[], double *vector, int n, double stopVal) {
    // tworzenie macierzy D^-1 i L+U, czyli macierzy iteracji iMatrix
    double ** dMatrix = new double*[n];
    double ** luMatrix = new double*[n];
    for (int i = 0; i < n; i++) { dMatrix[i] = new double[n]; luMatrix[i] = new double[n]; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            luMatrix[i][j] = matrix[i][j];
        }
    }
    for (int i = 0; i < n; i++) { dMatrix[i][i] = ((-1)/matrix[i][i]); luMatrix[i][i] = 0.0; }

    double **iterMatrix = multiplyMatrixByMatrix(n, dMatrix, luMatrix);
    for (int i = 0; i < n; i++) { delete[] dMatrix[i]; delete[] luMatrix[i]; }
    delete[] dMatrix; delete[] luMatrix;

    double *endVector, product;
    bool end = false;
    int iter = 0;
    while(!end) {
        endVector = multiplyMatrixByVector(n, n, (const double **)iterMatrix, vector);
        product = vectorEuclideanNorm(endVector, n);
        for (int i = 0; i < n; i++) { endVector[i] = endVector[i]/product; }

        if ( fabs(vectorEuclideanNorm(endVector, n) - vectorEuclideanNorm(vector, n)) < stopVal ) {
            end = true;
        } else {
            delete[] vector;
            vector = endVector;
            endVector = NULL;
        }
        iter++;
    }
    // cout << "Ilość iteracji: " << iter << endl;

    double retVal = fabs(rayleighQuotient(iterMatrix, endVector, n));

    for (int i = 0; i < n; i++) { delete[] iterMatrix[i]; }
    delete[] iterMatrix;

    return retVal;
}

double rayleighQuotient(double *matrix[], double *vector, int n) {
    double *tmpVector = multiplyMatrixByVector(n, n, (const double **)matrix, vector);

    double upper = 0, lower = 0;
    for (int i = 0; i < n; i++) {
        upper += (vector[i] * tmpVector[i]);
        lower += (vector[i] * vector[i]);
    }

    return upper/lower;
}

double *SORAlgorithm(int n, const double *matrix[], double *initVectorX, double *vectorB, double ro, double relax, int stopType) {
    bool finish = false;
    double *nextVector = new double[n];
    double sum, b, *tmpSwitch;
    int iter;
    for (iter = 1; !finish; iter++) {

        for (int i = 0; i < n; i++) {
            // sumowanie wiersz * wektor X
            sum = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    if (j < i) {
                        sum += matrix[i][j] * nextVector[j];
                    } else {
                        sum += matrix[i][j] * initVectorX[j];
                    }
                }
            }

            b = 0;
            b = vectorB[i]/matrix[i][i];
            sum = sum/matrix[i][i];
            nextVector[i] = b - sum;
        }

        // relaksacja
        for (int i = 0; i < n; i++) {
            nextVector[i] = ((1 - relax)*initVectorX[i]) + (relax * nextVector[i]);
        }

        if (stopType == 1) {
            finish = firstStopCriterium(ro, initVectorX, nextVector, n);
        } else {
            finish = secondStopCriterium(ro, nextVector, vectorB, matrix, n);
        }

        if (!finish) {
            tmpSwitch = nextVector;
            nextVector = initVectorX;
            initVectorX = tmpSwitch;
            for (int i = 0; i < n; i++) nextVector[i] = 0;
        }

    }

    cout << "Przeprowadzono " << iter << " iteracji" << endl;

    return nextVector;
}