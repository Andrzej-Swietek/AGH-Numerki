//
// Created by andrzej on 06.05.24.
//
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f(double x) {
    return log(x * x * x * x * x + 3 * x * x + x + 9);
}

double golden_ratio_minimization(double xa, double xb, double epsilon) {
    const double r = (sqrt(5.0) - 1) / 2;
    double x1, x2, x_min;
    int iteration = 0;
    ofstream outputFile("data/golden_ratio_output.txt");

    do {
        x1 = xa + r * (xb - xa);
        x2 = xa + (1 - r) * (xb - xa);

        if (f(x1) < f(x2)) {
            xb = x2;
        } else {
            xa = x1;
        }

        x_min = (x1 + x2) / 2;
        double diff = abs(-0.1673198 - x_min);

        outputFile << iteration << "\t" << x_min << "\t" << diff << endl;

        iteration++;
    } while (abs(x1 - x2) >= epsilon);

    outputFile.close();

    return x_min;
}

double trisection_minimization(double xa, double xb, double epsilon) {
    const double lambda1 = 1.0 / 3.0;
    const double lambda2 = 2.0 / 3.0;
    double x1, x2, x_min;
    int iteration = 0;
    ofstream outputFile("data/trisection_output.txt");

    do {
        x1 = xa + lambda1 * (xb - xa);
        x2 = xa + lambda2 * (xb - xa);

        if (f(x1) < f(x2)) {
            xb = x2;
        } else {
            xa = x1;
        }

        x_min = (x1 + x2) / 2;
        double diff = abs(-0.1673198 - x_min);

        outputFile << iteration << "\t" << x_min << "\t" << diff << endl;

        iteration++;
    } while (abs(x1 - x2) >= epsilon);

    outputFile.close();

    return x_min;
}

int main() {
    double xa = -0.5, xb = 1.0, epsilon = 1e-6;

    // Golden ratio minimization
    cout << "Golden Ratio Minimization:" << endl;
    double golden_min = golden_ratio_minimization(xa, xb, epsilon);
    cout << "Minimum found at x = " << golden_min << endl;

    // Trisection minimization
    cout << "\nTrisection Minimization:" << endl;
    double trisection_min = trisection_minimization(xa, xb, epsilon);
    cout << "Minimum found at x = " << trisection_min << endl;

    return 0;
}

