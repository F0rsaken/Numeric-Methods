#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include "../../lib/matrix-lib.h"
using namespace std;

const double xTab[] = {1.0, -1.0};
int n;
double **mainMatrix, *xVector, *answer, *freeColumn, *solution;
bool countedWithThomas;

void fillMatrix () {
    // k = 6, m = 4
    double k = 6.0, m = 4.0;
    mainMatrix = new double*[3];
    for (int i = 0; i < 3; i++) { mainMatrix[i] = new double[n]; }

    for (int i = 0; i < n; i++) { mainMatrix[1][i] = k; }
    for (int i = 0; i < n-1; i++) { mainMatrix[0][i] = 1.0/(i+1+m); }
    mainMatrix[0][n-1] = 0;
    mainMatrix[2][0] = 0;
    for (int i = 1; i < n; i++) { mainMatrix[2][i] = k/(i+2+m); }
}

void fillMatrixNormal () {
    // k = 6, m = 4
    double k = 6.0, m = 4.0;
    mainMatrix = new double*[n];
    for (int i = 0; i < n; i++) { mainMatrix[i] = new double[n]; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mainMatrix[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++) { mainMatrix[i][i] = k; }
    for (int i = 0; i < n-1; i++) { mainMatrix[i][i+1] = 1.0/(i+1+m); }
    for (int i = 1; i < n; i++) { mainMatrix[i][i-1] = k/(i+2+m); }
}

void clean() {
    if (countedWithThomas) {
        for (int i = 0; i < 3; i++) { delete[] mainMatrix[i]; }
    } else {
        for (int i = 0; i < n; i++) { delete[] mainMatrix[i]; }
    }
    delete[] mainMatrix;
    delete[] xVector;
    delete[] answer;
    delete[] solution;
}

int main( int argc, const char *argv[] ) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 0;
    }

    n = atoi(argv[1]);
    if (atoi(argv[2]) == 0) {
        countedWithThomas = false;
    } else {
        countedWithThomas = true;
    }
    srand(time(NULL));

    if (countedWithThomas) {
        cout << "Używany jest algorytm Thomasa" << endl;
        fillMatrix();
    } else {
        cout << "Używany jest algorytm Gaussa" << endl;
        fillMatrixNormal();
    }

    // cout << "Main matrix: \n";
    // printMatrix(3, n, mainMatrix);

    xVector = generateXVector(n, xTab, 2);

    auto start = chrono::system_clock::now(); 
    if (countedWithThomas) {
        answer = multiplyTridiagonalMatrixByVector(n, mainMatrix, xVector);
        solution = thomasAlgorithm(mainMatrix[2], mainMatrix[1], mainMatrix[0], answer, n);
    } else {
        answer = multiplyMatrixByVector(n, n, mainMatrix, xVector);
        solution = gaussElimination(n, n, mainMatrix, answer);
    }
    auto end = chrono::system_clock::now();
    // cout << "X vector: \n";
    // printVector(n, xVector);
    // cout << "Answer: \n";
    // printVector(n, answer);

    // cout << "Counted solution: \n";
    // printVector(n, solution);
    double xVectorLength = vectorEuclideanNorm(xVector, n);
    double solutionNorm = vectorEuclideanNorm(solution, n);
    chrono::duration<double> elapsed = end - start;

    cout << "Norma początkowego wektora X: " << xVectorLength << endl;
    cout << "Norma obliczonego wektora X: " << solutionNorm << endl;
    cout << fixed << setprecision(3) << "Różnica procentowa wektorów: " << abs(1 - (solutionNorm / xVectorLength)) << '%' << endl;

    cout << "Czas trwania algorytmu: " << elapsed.count() << "s\n";

    clean();

    return 0;
}