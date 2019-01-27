#include <iostream>
#include <fstream>
#include <math.h>
#include "matrix-lib.h"
#include "diff-equations-lib.h"
#include "interpolation-lib.h"
using namespace std;

void sendPlotToFileDiff(PointDifferential data[], int n, string fileName, bool informUser) {
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

Point * MRSAlgorithm(double h, int n, Point startP, Point endP) {
    double * upper = new double[n-1];
    double * lower = new double[n-1];
    double * a = new double[n];
    double * d = new double[n];
    Point *outPoints = new Point[n];

    double factor = 1/pow(h, 2);

    for (int i = 0; i < n-1; i++) {
        upper[i] = factor;
        lower[i] = factor;
    }
    upper[0] = 0;
    lower[n-2] = 0;

    outPoints[0].x = startP.x;
    outPoints[n - 1].x = endP.x;
    for (int i = 1; i < n-1; i++) { outPoints[i].x = outPoints[i-1].x + h; }

    d[0] = startP.y;
    d[n-1] = endP.y;
    for (int i = 1; i < n-1; i++) {
        a[i] = 4 - (2*factor);
        d[i] = 4*outPoints[i].x;
    }

    a[0] = 1;
    a[n-1] = 1;

    double *y = thomasAlgorithm(lower, a, upper, d, n);

    for (int i = 0; i < n; i++) {
        outPoints[i].y = y[i];
    }

    delete[] upper;
    delete[] lower;
    delete[] a;
    delete[] d;

    return outPoints;
}