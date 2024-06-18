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
    if (n % 3 != 0) {
        throw invalid_argument("n must be a multiple of 3 for Simpson's 3/8 rule.");
    }
    double sum = f[0] + f[n];
    for (int i = 1; i < n; i++) {
        if (i % 3 == 0) {
            sum += 2 * f[i];
        } else {
            sum += 3 * f[i];
        }
    }
    return (3 * h / 8) * sum;
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

vector<vector<double>> richardsonExtrapolation(const function<double(double)>& f, double a, double b, int maxN, bool useSimpson) {
    vector<vector<double>> integrals(maxN + 1, vector<double>(maxN + 1));

    for (int n = 0; n <= maxN; n++) {
        int N = useSimpson ? 3 * pow(2, n) : pow(2, n);
        double h = (b - a) / N;
        vector<double> f_values(N + 1);

        for (int i = 0; i <= N; ++i) {
            double x = a + i * h;
            f_values[i] = f(x);
        }

        if (useSimpson) {
            integrals[n][0] = simpsonIntegral(f_values, h, N);
        } else {
            integrals[n][0] = trapezoidalIntegral(f_values, h, N);
        }

        for (int k = 1; k <= n; k++) {
            integrals[n][k] = (pow(4, k) * integrals[n][k - 1] - integrals[n - 1][k - 1]) / (pow(4, k) - 1);
        }
    }

    return integrals;
}


void saveResults(const vector<vector<double>>& integrals, const string& filename) {
    ofstream file(filename);
    for (const auto& row : integrals) {
        for (const auto& value : row) {
            file << value << " ";
        }
        file << endl;
    }
    file.close();
}


int main() {
    double a = 0;
    double b = 1;
    int maxN = 8;

    // Define the function f(x)
    function<double(double)> func = [](double x) {
        return log(x * x * x + 3 * x * x + x + 0.1) * sin(18 * x);
    };

    // Calculate integrals using Richardson Extrapolation
    auto integralsTrapezoidal = richardsonExtrapolation(func, a, b, maxN, false);
    auto integralsSimpson = richardsonExtrapolation(func, a, b, maxN, true);

    // Save results to files
    saveResults(integralsTrapezoidal, "data/integrals_trapezoidal.txt");
    saveResults(integralsSimpson, "data/integrals_simpson.txt");

    // Print results for Trapezoidal rule
    cout << "Trapezoidal Rule:" << endl;
    for (int n = 0; n <= maxN; n++) {
        cout << "N = " << n << ": " << integralsTrapezoidal[n][n] << endl;
    }

    // Print results for Simpson's 3/8 rule
    cout << endl << "Simpson's 3/8 Rule:" << endl;
    for (int n = 0; n <= maxN; n++) {
        cout << "N = " << n << ": " << integralsSimpson[n][n] << endl;
    }


    return 0;
}
