#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "../../lib/matrix-lib.h"
using namespace std;

const double xTab[] = {1.0, -1.0};
int n;
double **mainMatrix, * xVector, * answer, * freeColumn, * solution;

void fillMatrix () {
    mainMatrix = new double*[n];
    int i;
    for (i = 0; i < n; i++) { mainMatrix[i] = new double[n]; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                mainMatrix[i][j] = (2*(i+1))/(j+1);
            } else {
                mainMatrix[i][j] = mainMatrix[j][i];
            }
        }
    }

    // for (i = 0; i < n; i++) { mainMatrix[0][i] = 1.0; }

    // for (i = 1; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
    //         mainMatrix[i][j] = 1.0 / (i + j + 1.0);
    //     }
    // }
}

void clean() {
    for (int i = 0; i < n; i++) { delete[] mainMatrix[i]; }
    delete[] mainMatrix;
    delete[] xVector;
    delete[] answer;
    delete[] solution;
}

int main( int argc, const char *argv[] ) {
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 0;
    }

    n = atoi(argv[1]);
    srand(time(NULL));

    fillMatrix();
    // cout << "Main matrix: \n";
    // printMatrix(n, n, mainMatrix);

    xVector = generateXVector(n, xTab, 2);
    answer = multiplyMatrixByVector(n, n, mainMatrix, xVector);
    // cout << "X vector: \n";
    // printVector(n, xVector);
    // cout << "Answer: \n";
    // printVector(n, answer);

    solution = gaussElimination(n, n, mainMatrix, answer);
    // cout << "Counted solution: \n";
    // printVector(n, solution);
    double xVectorLength = vectorEuclideanNorm(xVector, n);
    double solutionNorm = vectorEuclideanNorm(solution, n);

    cout << "Norma początkowego wektora X: " << xVectorLength << endl;
    cout << "Norma obliczonego wektora X: " << solutionNorm << endl;
    cout << fixed << setprecision(3) << "Różnica procentowa wektorów: " << abs(1 - (solutionNorm/xVectorLength)) << '%' << endl;

    clean();

    return 0;
}