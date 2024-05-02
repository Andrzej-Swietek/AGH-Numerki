#include<iostream>
#include<fstream>
#include<sstream>
#include<math.h>
#include<vector>
#include<functional>

#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"
#include "/usr/include/gsl/gsl_eigen.h"

using namespace std;

void calculate_derivatives(double* derivatives, int n){
    for (int k = 0; k <= n; ++k) {
        if(k % 2 == 0) {
            derivatives[k] = pow(-1, k/2);
        } else {
            derivatives[k] = 0;
        }
    }
}

constexpr long long factorial(int n) {
    return n <= 1 ? 1 : n * factorial(n-1);
}

// Funkcja rozwiązująca układ równań
void solve_system(double* c, double* f_k, double* b, int n) {
    gsl_matrix* A = gsl_matrix_alloc(n, n);
    gsl_vector* x = gsl_vector_alloc(n);
    gsl_vector* y = gsl_vector_alloc(n);
    gsl_permutation* p = gsl_permutation_alloc(n);

    // Inicjalizacja macierzy A i wektora y
    for (int i = 0; i < n; ++i) {
        gsl_vector_set(y, i, -f_k[i]);
        for (int j = 0; j < n; ++j) {
            gsl_matrix_set(A, i, j, c[i+j]);
        }
    }
    
    // Rozwiązanie układu równań
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, y, x);

    // Przypisanie rozwiązania do wektora b
    for (int i = 0; i < n; ++i) {
        b[i] = gsl_vector_get(x, i);
    }

    // Zwolnienie pamięci
    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_permutation_free(p);
}

double R_NM(double x, double *a, double *b, int N, int M){
    double sum_a = 0.0;
    double sum_b = 0.0;
    for (int i = 0; i <= N; ++i) sum_a += a[i] * pow(x, i);
    for (int i = 0; i <= M; ++i) sum_b += b[i] * pow(x, i);
    return sum_a / sum_b;
}

int main() {
    // Stopień pochodnych
   int n = 10;

    // Inicjalizacja tablicy na przechowywanie pochodnych
    double f_k[n+1];

    calculate_derivatives(f_k, n);

    cout << "Pochodne funkcji cos(x) w punkcie x = 0:\n";
    for (int k = 0; k <= n; ++k) {
        cout << "f^(" << k << ")(0) = " << f_k[k] << endl;
    }

    cout << "\nc_k:\n";
    for (int k = 0; k <= n; ++k) {
        f_k[k] = f_k[k] / factorial(k); // od tego momentu f_k to c_k
        cout << "c_" << k << " = " << f_k[k] << endl;
    }

    // Stopnie wielomianów Q_M(x) i P_N(x)
    int M = 4;
    int N = 4;

    // Inicjalizacja tablic na współczynniki wielomianów
    double c[N+M+1], b[M+1];

    // Rozwiązanie układu równań
    solve_system(c, f_k, b, N+M+1);

    // Wyświetlanie współczynników wielomianów
    cout << "\nWspółczynniki wielomianu Q_M(x):\n";
    for (int i = 0; i <= M; ++i) 
        cout << "b_" << i << " = " << b[i] << endl;

    
    ofstream file("data/interpolation_results_" + to_string(N) + "_" + to_string(M) + ".csv");
    for (double x = -5.0; x <= 5.0; x += 0.1 ){
        double interpolated = R_NM(x, f_k, b, N, M);
        file << x << " " << interpolated << endl;
    }
    file.close();

    return 0;
}
