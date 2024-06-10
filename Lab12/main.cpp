#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

using namespace std;
using dbl = double;

double f(double x) { return log( pow(x,3) + 3*pow(x,2) + x + 0.1 )*sin(18*x); }


// double trapezoid_integration(vector<double>f, int N, double h ) {
//     double sum =0.0;
//     for( int i = 0; i < N-1; i++ ) {
//         sum += 0.5*h *( f[i]  + f[i+1]); 
//     }

//     return sum;
// }
// double three_eighth(vector<double>f, int N, double h){
//     double sum = 0.0;
//      for( int i = 0; i < (N/3)-1; i++ ) 
//         sum += 3*h/8 * ( f[3*i] + 3* f[3*i+1] + 3*f[3*i+2] + f[3*i+3] );

//     return sum;
// }



// Function to calculate the Trapezoidal rule integral
double trapezoidalIntegral(const vector<double>& f, double h, int n) {
    double sum =0.0;
        for( int i = 0; i <= n+1; i++ ) {
        sum += 0.5*h *( f[i]  + f[i+1]); 
    }
    return sum;
}

// Function to calculate the Simpson's 3/8 rule integral
double simpsonIntegral(const vector<double>& f, double h, int n) {
    double sum1 = 0, sum2 = 0;
    for (int i = 1; i < n; i++) {
        sum1 += f[i];
        sum2 += f[2 * i];
    }
    return h * (f[0] + f[n] + 4 * sum1 + 2 * sum2) / 3;
}

// Function to perform Richardson Extrapolation
vector<vector<double>> richardsonExtrapolation(const vector<double>& f, double a, double b, int maxN) {
    vector<vector<double>> integrals(maxN + 1, vector<double>(maxN + 1));

    // Calculate integrals using Trapezoidal rule for different step sizes
    for (int n = 1; n <= maxN; n++) {
        double h = (b - a) / n;
        integrals[n][0] = trapezoidalIntegral(f, h, n);

        for (int k = 1; k <= n; k++) {
            // Perform Richardson Extrapolation
            integrals[n][k] = (4 * integrals[n][k - 1] - integrals[n - 1][k - 1]) / (4 - pow(2, -k));
        }
    }

    return integrals;
}

int main() {
    // Define the function f(x)
    function<double(double)> f = [](double x) {
        return log(x * x * x + 3 * x * x + x + 0.1) * sin(18 * x);
    };

    // Define the integration limits
    double a = 0;
    double b = 1;

    // Set the maximum number of nodes
    int maxN = 8;

    // Define discretization (number of points)
    int N = 100;

    // Create a vector to store function values
    vector<double> f_values(N + 1);

    // Calculate function values at each point
    for (int i = 0; i <= N; ++i) {
        double x = a + (b - a) * i / N;
        f_values[i] = f(x);
    }

    // Call richardsonExtrapolation with the pre-computed vector
    vector<vector<double>> integralsTrapezoidal = richardsonExtrapolation(f_values, a, b, maxN);
    vector<vector<double>> integralsSimpson = richardsonExtrapolation(f_values, a, b, maxN);

    // Print the results for Trapezoidal rule
    cout << "Trapezoidal Rule:" << endl;
    for (int n = 0; n <= maxN; n++) {
        cout << "N = " << n << ": " << integralsTrapezoidal[n][n] << endl;
    }

    // Print the results for Simpson's 3/8 rule
    cout << endl << "Simpson's 3/8 Rule:" << endl;
    for (int n = 0; n <= maxN; n++) {
        cout << "N = " << n << ": " << integralsSimpson[n][n] << endl;
    }

    for ( auto row :  integralsTrapezoidal) {
        for ( auto column : row ) 
            cout << column << " ";
        cout << "\n";
    }

    for ( auto row :  integralsSimpson) {
        for ( auto column : row ) 
            cout << column << " ";
        cout << "\n";
    }


    return 0;
}
