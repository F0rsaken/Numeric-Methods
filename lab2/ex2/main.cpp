#include <cstdlib>
#include <iostream>
#include <time.h>
#include "../../lib/matrix-lib.h"
using namespace std;

const double xTab[] = {1.0, -1.0};
double **mainMatrix;
int n;
double ro = 0.01;

// algorytm jest zbieżny dla promienia spektralnego macierzy D^-1(L + U) < 1

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
    // delete[] countedVector;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    n = atoi(argv[1]);
    srand(time(NULL));

    fillMatrix();
    // printMatrix(n, n, mainMatrix);
    double * initVector = new double[n];
    for (int i = 0; i < n; i++) initVector[i] = 1;
    double radius = powerIteration(mainMatrix, initVector, n, ro);
    cout << "Promień spektralny macierzy: " << radius << endl;

    clean();
    delete[] initVector;

    return 0;
}