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

void rungeKuttyMethodDiff(PointDifferential points[], double h, int n, double yPrim(double x, double y)) {
    double k1, k2, k3, k4, dyn;

    for (int i = 1; i < n; i++) {
        points[i].x = points[i-1].x + h;
        k1 = h*points[i-1].dY;
        k2 = h*yPrim(points[i-1].x + (h/2), points[i-1].y + (k1/2));
        k3 = h*yPrim(points[i-1].x + (h/2), points[i-1].y + (k2/2));
        k4 = h*yPrim(points[i-1].x + h, points[i-1].y + k3);
        dyn = (k1 + (2*k2) + (2*k3) + k4)/6;
        points[i].y = points[i-1].y + dyn;
        points[i].dY = yPrim(points[i].x, points[i].y);
    }
}