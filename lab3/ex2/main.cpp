#include <cstdlib>
#include <math.h>
#include <iostream>
#include "../../lib/equations-lib.h"
#include "../../lib/matrix-lib.h"
using namespace std;

bool stop0 = true;
double ro = 0.1;


int main(int argc, char const *argv[])
{
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 1;
        // cout << "Błąd! Nie potrzeba argumentów\nWszystkie argumenty będą zignorowane" << endl;
    }

    int stopCriterium = atoi(argv[1]);

    FunctionOf3 *fVector = new FunctionOf3[3];
    fVector[0].setF( [](double x1, double x2, double x3) { return pow(x1, 2) + pow(x2, 2) + x3 - 1; } );
    fVector[1].setF( [](double x1, double x2, double x3) { return (2 * pow(x1, 2)) - pow(x2, 2) - (4 * pow(x3, 2)) + 3; } );
    fVector[2].setF( [](double x1, double x2, double x3) { return pow(x1, 2) + x2 + x3 - 1; } );

    FunctionOf1 **jacobianMatrix = new FunctionOf1*[3];
    for (int i = 0; i < 3; i++) { jacobianMatrix[i] = new FunctionOf1[3]; }

    jacobianMatrix[0][0].setF( [](double x) { return 2*x; } );
    jacobianMatrix[0][1].setF( [](double x) { return 2*x; } );
    jacobianMatrix[0][2].setF( [](double x) { return 1.0; } );

    jacobianMatrix[1][0].setF( [](double x) { return 4 * x; } );
    jacobianMatrix[1][1].setF( [](double x) { return -2 * x; } );
    jacobianMatrix[1][2].setF( [](double x) { return -8 * x; } );

    jacobianMatrix[2][0].setF( [](double x) { return 2 * x; } );
    jacobianMatrix[2][1].setF( [](double x) { return 1.0; } );
    jacobianMatrix[2][2].setF( [](double x) { return 1.0; } );

    double initVector[3] = {1, 1, 1};
    double * solution = newtonEquationsSystem(fVector, jacobianMatrix, initVector, stopCriterium, ro);
    cout << "Rozwiązanie x1 = " << solution[0] << ", x2 = " << solution[1] << ", x3 = " << solution[2] << endl;

    delete[] solution;
    for (int i = 0; i < 3; i++) { delete[] jacobianMatrix[i]; }
    delete[] jacobianMatrix;


    return 0;
}