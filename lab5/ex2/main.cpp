#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>
#include "../../lib/matrix-lib.h"
#include "../../lib/interpolation-lib.h"
using namespace std;

double step = 0.1;

Point * getDataPoints(int n) {
    
}

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        cout << "Zła liczba argumentów!\n";
        return 1;
    }

    // int printMode = atoi(argv[1]); // 0 - cout, 1 - do pliku
    // int boundaryType = atoi(argv[1]); // 0 - lagrange, 1 - Newton
    int n = atoi(argv[1]);
    
    // if (boundaryType != 0 && boundaryType != 1) {
        // cout << "Zła wartość pierwszego argumentu!\n";
        // return 1;
    // }

    Point *data = getDataPoints(n);
    Polynomial2 *splines = quadraticSplines(data, n);
    
    // int outN = ((data[n - 1].x - data[0].x) / step) + 1;
    // cout << "OutN: " << outN << endl;
    // Point *interpolation = new Point[outN];
    
    // cout << "Interpolacja: \n";
    // ofstream outputFile;
    // outputFile.open("out.dat", ios::trunc);
    // for(int i = 0; i < outN; i++) {
    //     outputFile << interpolation[i].x << " " << interpolation[i].y << endl;
    //     // cout << "x: " << interpolation[i].x << ", y: " << interpolation[i].y << endl;
    // }
    // outputFile.close();
    // cout << "Skończono pisać do pliku\n";

    delete[] data;
    delete[] splines;

    return 0;
}
