#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <functional>
#include <iomanip>

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

void print_gsl_matrix(gsl_matrix *A) {
    if (A == nullptr) {
        std::cout << "Matrix is empty." << std::endl;
        return;
    }

    int rows = A->size1;
    int cols = A->size2;

    std::cout << "Matrix[" << rows << "][" << cols << "]:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::setw(10) << gsl_matrix_get(A, i, j) << " ";
        }
        std::cout << std::endl;
    }
}


// Funkcja rozwiązująca układ równań
void solve_system(double* c, double* f_k, double* b, int n, int N, int M) {
    gsl_matrix* A = gsl_matrix_alloc(M, M);
    gsl_vector* x = gsl_vector_alloc(M);
    gsl_vector* y = gsl_vector_alloc(M);
    gsl_permutation* p = gsl_permutation_alloc(M);

    // Inicjalizacja macierzy A i wektora y
    for (int i = 0; i < M; ++i) {
        gsl_vector_set(y, i, -f_k[N - M + i + 1]);
        for (int j = 0; j < M; ++j) {
            gsl_matrix_set(A, i, j, c[N - M + i + j + 1]);
        }
    }

    print_gsl_matrix(A);

    // Rozwiązanie układu równań
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, y, x);

    for( int i = 0; i < M; i++ )
        b[i] = gsl_vector_get(x, i);

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
    return sum_b == 0? 0 :  sum_a / sum_b;
}

void calculate_vector_c(double*c, double* f_k, int n) {
    cout << "\nc_k:\n";
    for (int k = 0; k <= n; ++k) {
        c[k] = f_k[k] / factorial(k); // od tego momentu f_k to c_k
        cout << "c_" << k << " = " << c[k] << endl;
    }
}

void calculate_vector_a(double*a, double* b, double* c, int N, int M) {
    fill(a, a + N + 1, 0.0);
    for (int i = 0; i <= N; i++ )
        for (int j =0; j <= i; j++)
            a[i] += c[i-j] * b[j];
}

int main() {

    // Stopnie wielomianów Q_M(x) i P_N(x)
    int M = 4;
    int N = 4;

    // Stopień pochodnych
    int n = M+N;

    // Inicjalizacja tablicy na przechowywanie pochodnych
    double f_k[n+1];
    double c[n+1];

    calculate_derivatives(f_k, n);

    cout << "Pochodne funkcji cos(x) w punkcie x = 0:\n";
    for (int k = 0; k <= n; ++k) {
        cout << "f^(" << k << ")(0) = " << f_k[k] << endl;
    }

    calculate_vector_c(c, f_k, n);


    // Inicjalizacja tablic na współczynniki wielomianów
    double matrix_c[N+M+1], b[M+1];
    for (int i = 0; i <= M-1; i++)
        for(int j = 0; j <= M-1; j++ )
            matrix_c[ i + j ] = c[N - M + i +1 ]; // lub   matrix_c[i + j] = c[j]; etc

    // Rozwiązanie układu równań
    solve_system(matrix_c, c, b, n, N, M);

    // Wyświetlanie współczynników wielomianów
    cout << "\nWspółczynniki wielomianu Q_M(x):\n";
    for (int i = 0; i <= M; ++i)
        cout << "b_" << i << " = " << b[i] << endl;

    // vector a
    double a[n+1];
    calculate_vector_a(a, b, c, N, M);

    cout << "\nWspółczynniki wielomianu P_N(x):\n";
    for (int i = 0; i <= N; ++i)
        cout << "a_" << i << " = " << a[i] << endl;

    ofstream file("data/approximation_results_" + to_string(N) + "_" + to_string(M) + ".csv");
    file << "x y\n";
    for (double x = -5.0; x <= 5.0; x += 0.1 ){
        double approximated = R_NM(x, a, b, N, M);
        file << x << " " << approximated << endl;
    }
    file.close();

    return 0;
}
