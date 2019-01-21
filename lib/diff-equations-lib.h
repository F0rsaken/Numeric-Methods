#ifndef DIFF_EQUATIONS_LIB
#define DIFF_EQUATIONS_LIB
#include <string>
using namespace std;

// class Point {
// public:
//     double x;
//     double y;

//     Point() {}
//     Point(double x, double y) {
//         this->x = x;
//         this->y = y;
//     }
// };

class PointDifferential {
public:
    double x;
    double y;
    double dY;

    PointDifferential() {}
    PointDifferential(double x, double y, double dY) {
        this->x = x;
        this->y = y;
        this->dY = dY;
    }
    ~PointDifferential() {}
};

void sendPlotToFileDiff(PointDifferential data[], int n, string fileName, bool informUser = false);

void eulerMethodDiff(PointDifferential points[], double h, int n, double yPrim(double x, double y));
void rungeKuttyMethodDiff(PointDifferential points[], double h, int n, double yPrim(double x, double y));

#endif