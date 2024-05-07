#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;


double f_1(double x) {
    return log( pow(x,5) + 3* pow(x,2) + x + 9 );
}


double f_2(double x) {
    return pow(x,6);
}


double F(double x1, double x2, double (*f)(double)) {
    return (f(x2) - f(x1)) / (x2 - x1);
}


double F2(double x1, double x2, double x3, double (*f)(double)) {
    return  (
        (f(x3) - f(x2)) / (x3 - x2) - (f(x2) - f(x1)) / (x2 - x1) 
     ) / (x3 - x1);
}


double powell_interpolation(double x1, double x2, double x3, double (*f)(double)) {
    return (x1 + x2)/2 - 0.5 * F(x1, x2, f) / F2(x1, x2, x3, f);
}


void find_minimum(double (*f)(double), double x1, double x2, double x3, int iterations, double h, const std::string& filename) {
    std::ofstream file(filename);
    file << "x1\tx2\tx3\txm\tF[x1,x2]\tF[x1,x2,x3]\n";
    for (int i = 0; i < iterations; ++i) {
        double xm = powell_interpolation(x1, x2, x3, f);
        file << x1 << "\t" << x2 << "\t" << x3 << "\t" << xm << "\t" << F(x1, x2, f) << "\t" << F2(x1, x2, x3, f) << "\n";
 
        double max = fabs(x1-xm);
        if(fabs(x2-xm) > max){
            max = fabs(x2-xm);
            x2 = xm;
        }
        else if(fabs(x3-xm) > max){
            max = fabs(x3-xm);
            x3 = xm;
        }
        else{
            x1 = xm;
        }
    }
}


int main() {
    // double epsilon = 1e6;
    double h = 0.01;

    // double x_min = -1.5; double x_max = 1.0;

    double x_1 = -0.5;
    double x_2 = x_1 + h;
    double x_3 = x_2 + h; 
    
    find_minimum(f_1, x_1, x_2, x_3, 10, h, "data/f1_results.csv");

    //------------------------------------------------------------------------

    x_1 = -0.9;
    find_minimum(f_1, x_1, x_1 + h, x_1 + 2 * h, 10, h, "data/f1_results_2.csv");

    //------------------------------------------------------------------------

    double x_1_f2 = 1.5;
    double x_2_f2 = x_1 + h;
    double x_3_f2 = x_2 + h; 

    find_minimum(f_2, x_1_f2, x_2_f2, x_3_f2, 100, h, "data/f2_results.csv");

    return 0;
}