#include <iostream>
#include <fstream>
#include "diff-equations-lib.h"
using namespace std;

void sendPlotToFile(PointDifferential data[], int n, string fileName, bool informUser) {
    // cout << "Interpolacja: \n";
    ofstream outputFile;
    outputFile.open(fileName, ios::trunc);
    for(int i = 0; i < n; i++) {
        outputFile << data[i].x << " " << data[i].y << endl;
    }
    outputFile.close();
    if (informUser) cout << "Skończono pisać do pliku " << fileName << endl;
}

// pierwszy wiersz points musi byc usupelniony
void eulerMethodDiff(PointDifferential points[], double h, int n, double yPrim(double x, double y)) {
    for (int i = 1; i < n; i++) {
        points[i].x = points[i-1].x + h;
        points[i].y = points[i-1].y + (h * points[i-1].dY);
        points[i].dY = yPrim(points[i].x, points[i].y);
    }
}