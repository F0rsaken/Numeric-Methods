#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "../../lib/matrix-lib.h"
using namespace std;

const float xTab[] = {1.0, -1.0};
int n;
float **mainMatrix, * xVector, * answer, * freeColumn, * solution;

void fillMatrix () {
    mainMatrix = new float*[n];
    int i;
    for (i = 0; i < n; i++) { mainMatrix[i] = new float[n]; }

    for (i = 0; i < n; i++) { mainMatrix[0][i] = 1.0; }

    for (i = 1; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mainMatrix[i][j] = 1.0 / (i + j + 1.0);
        }
    }
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

    xVector = generateXVector_f(n, xTab, 2);
    answer = multiplyMatrixByVector_f(n, n, mainMatrix, xVector);
    // cout << "X vector: \n";
    // printVector(n, xVector);
    // cout << "Answer: \n";
    // printVector(n, answer);

    solution = gaussElimination_f(n, n, mainMatrix, answer);
    // cout << "Counted solution: \n";
    // printVector(n, solution);
    float xVectorLength = vectorEuclideanNorm_f(xVector, n);
    float solutionNorm = vectorEuclideanNorm_f(solution, n);

    cout << "Norma początkowego wektora X: " << xVectorLength << endl;
    cout << "Norma obliczonego wektora X: " << solutionNorm << endl;
    cout << fixed << setprecision(3) << "Różnica procentowa wektorów: " << abs(1 - (solutionNorm/xVectorLength)) << '%' << endl;

    clean();

    return 0;
}