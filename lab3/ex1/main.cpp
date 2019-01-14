#include <cstdlib>
#include <math.h>
#include <iostream>
#include "../../lib/equations-lib.h"
// #include "../../lib/matrix-lib.h"
using namespace std;

bool stop0 = true;
double ro = 0.000001;
// f(x) = x^n - (1 - x)^m, gdzie n = 12, m = 15, a [a, b] = [-1.8, 0.6]
// f(x) = x^12 - (1 -x)^15
// f'(x) = n*x^(n-1) + m*(1-x)^(m-1)
// f'(x) = 12 * x^11 + 15 * (1-x)^14
// ro = (okolo) 0.1
// realRoot = 0.461403
// xi+1 = xi - F(xi)/F'(xi)

double mainFunction(double x) {
    double a = pow(x, 12);
    double b = pow(1-x, 15);
    return a - b;
}

double functionPrim(double x) {
    double a = 12*pow(x, 11);
    double b = 15*pow(1-x, 14);
    return a + b;
}

// double useNewton(double x0, int iteration) {
//     double fXi = mainFunction(x0);
//     double fpXi = functionPrim(x0);
//     double xK = x0 - (fXi/fpXi);
//     if (stop0) {
//         if ( fabs(xK-x0) < ro ) {
//             cout << "Iterations: " << iteration << endl;
//             return xK;
//         } else {
//             return useNewton(xK, iteration+1);
//         }
//     } else {
//         if ( mainFunction(xK) < ro ) {
//             cout << "Iterations: " << iteration << endl;
//             return xK;
//         } else {
//             return useNewton(xK, iteration+1);
//         }
//     }
// }

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Zła liczba argumentów!\n";
        return 1;
        // cout << "Błąd! Nie potrzeba argumentów\nWszystkie argumenty będą zignorowane" << endl;
    }

    // Argumenty: 1 | 2 - kryterium stopu, 1 | 2 - metoda Newtona | metoda Eulera

    int stopCriterium = atoi(argv[1]);
    int methodType = atoi(argv[2]);
    double a = -1.8;
    double b = 0.6;
    double rootVal;
    cout << "ro = " << ro << endl;
    if (methodType == 1) {
        cout << "Newton method: \n";
        for (double tmp = a; tmp < b; tmp += 0.1) {
            rootVal = newtonMethod(tmp, mainFunction, functionPrim, stopCriterium, ro);
            cout << "Wartość pierwiastka: " << rootVal << " dla startowego x = " << tmp << endl << endl;
        }
    } else if (methodType == 2) {
        cout << "Euler method: \n";
        for (double tmp = a; tmp < b; tmp += 0.1) {
            rootVal = eulerMethod(tmp, b, mainFunction, stopCriterium, ro);
            cout << "Wartość pierwiastka: " << rootVal << " dla startowego przedziału [a, b] = [" << tmp << ", " << b << "]\n" << endl;
        }
    } else {
        cout << "Zły drugi argument!\n";
    }
    return 0;
}