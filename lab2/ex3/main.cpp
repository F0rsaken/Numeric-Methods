#include <cstdlib>
#include <iostream>
#include <time.h>
#include "../../lib/matrix-lib.h"
using namespace std;

const double xTab[] = {1.0, -1.0};
double **mainMatrix, *vectorX, *vectorB, *countedSolution;
double ro = 0.01;
double omega = 1.5;
int n;

void fillMatrix() {
    mainMatrix = new double*[n];
    for (int i = 0; i < n; i++) { mainMatrix[i] = new double[n]; }

    double k = 5.0, m = 0.5;
    for (int i = 0; i < n; i++) { mainMatrix[i][i] = k; }

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            mainMatrix[i][j] = m/(j+1);
            if ( (j+1)%2 == 1 ) mainMatrix[i][j] *= -1;
        }
    }

    for (int i = 1; i < n; i++) { mainMatrix[i][i-1] = m/(i+1); }

    for (int i = 0; i < n; i++) {
        for (int j = i-2; j >= 0; j--) {
            mainMatrix[i][j] = 0;
        }
    }

}

void clean() {
    for (int i = 0; i < n; i++) {
        delete[] mainMatrix[i];
    }
    delete[] mainMatrix;
    delete[] vectorX;
    delete[] vectorB;
    delete[] countedSolution;
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    n = atoi(argv[1]);
    int stopTestType = atoi(argv[2]);
    if (stopTestType != 1 && stopTestType != 2) {
        cout << "Nieprawidłowy drugi argument! Dopuszczalne wartości: 1 | 2\n";
        return 1;
    }
    srand(time(NULL));

    fillMatrix();
    // printMatrix(n, n, mainMatrix);
    vectorX = generateXVector(n, xTab, 2);
    vectorB = multiplyMatrixByVector(n, n, (const double**)mainMatrix, vectorX);
    double * initVector = new double[n];
    for (int i = 0; i < n; i++) initVector[i] = 0;

    if (stopTestType == 1) {
        countedSolution = SORAlgorithm(n, (const double**) mainMatrix, initVector, vectorB, ro, omega, stopTestType);
    } else {
        countedSolution = SORAlgorithm(n, (const double **) mainMatrix, initVector, vectorB, ro, omega, stopTestType);
    }
    // printVector(n, vectorX);
    // printVector(n, countedSolution);
    double vectorXNorm = vectorEuclideanNorm(vectorX, n);
    double solutionNorm = vectorEuclideanNorm(countedSolution, n);

    cout << "Norma początkowego wektora X: " << vectorXNorm << endl;
    cout << "Norma obliczonego wektora X: " << solutionNorm << endl;

    clean();

    return 0;
}