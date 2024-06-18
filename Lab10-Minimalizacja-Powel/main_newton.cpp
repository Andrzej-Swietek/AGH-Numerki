//
// Created by andrzej on 06.05.24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

// Definicja funkcji g(x, y)
double g(double x, double y) {
    return x * x - 4 * x + y * y - 4 * y + x * y;
}

// Definicja pochodnych cząstkowych funkcji g(x, y)
double dg_dx(double x, double y) {
    return 2 * x - 4 + y;
}

double dg_dy(double x, double y) {
    return 2 * y - 4 + x;
}

// Definicja macierzy Hessego funkcji g(x, y)
gsl_matrix* hessian(double x, double y) {
    gsl_matrix* H = gsl_matrix_alloc(2, 2);
    gsl_matrix_set(H, 0, 0, 2);
    gsl_matrix_set(H, 0, 1, 1);
    gsl_matrix_set(H, 1, 0, 1);
    gsl_matrix_set(H, 1, 1, 2);
    return H;
}

// Definicja wektora gradientu funkcji g(x, y)
gsl_vector* gradient(double x, double y) {
    gsl_vector* grad = gsl_vector_alloc(2);
    gsl_vector_set(grad, 0, dg_dx(x, y));
    gsl_vector_set(grad, 1, dg_dy(x, y));
    return grad;
}

// Metoda Newtona do znalezienia minimum funkcji g(x, y)
pair<double, double> newton_method(double x, double y, double epsilon) {
    gsl_vector* r = gsl_vector_alloc(2);
    gsl_vector_set(r, 0, x);
    gsl_vector_set(r, 1, y);
    gsl_vector* delta_r = gsl_vector_alloc(2);
    gsl_matrix* H_inv = gsl_matrix_alloc(2, 2);
    gsl_permutation* p = gsl_permutation_alloc(2);
    int signum;
    int iterations = 0;
    double omega = 1.0; // Wartość domyślna wagi

    do {
        // Obliczenia macierzy odwrotnej Hessianu
        gsl_matrix* H = hessian(gsl_vector_get(r, 0), gsl_vector_get(r, 1));
        gsl_linalg_LU_decomp(H, p, &signum);
        gsl_linalg_LU_invert(H, p, H_inv);

        // Przepis iteracyjny
        gsl_blas_dgemv(CblasNoTrans, 1.0, H_inv, gradient(gsl_vector_get(r, 0), gsl_vector_get(r, 1)), 0.0, delta_r);
        gsl_vector_scale(delta_r, omega);
        gsl_vector_sub(r, delta_r);

        iterations++;
    } while (gsl_blas_dnrm2(delta_r) >= epsilon);

    double min_x = gsl_vector_get(r, 0);
    double min_y = gsl_vector_get(r, 1);

    gsl_vector_free(r);
    gsl_vector_free(delta_r);
    gsl_matrix_free(H_inv);
    gsl_permutation_free(p);

    return make_pair(min_x, min_y);
}

int main() {
    double epsilon = 1e-6;

    // Obliczanie minimum funkcji g(x, y)
    pair<double, double> min_point = newton_method(0, 0, epsilon);

    // Wypisywanie minimum funkcji
    cout << "Minimum of g(x, y) found at: x = " << min_point.first << ", y = " << min_point.second << endl;


    // Otwieranie pliku do zapisu danych
    ofstream outputFile("data/function_values_netwon_metd.txt");
    if (!outputFile.is_open()) {
        cerr << "Error: Unable to open output file." << endl;
        return 1;
    }

    // Obliczanie wartości funkcji g(x, y) w zakresie (-10, 10)
    for (double x = -10; x <= 10; x += 0.5) {
        for (double y = -10; y <= 10; y += 0.5) {
            double value = g(x, y);
            outputFile << x << " " << y << " " << value << endl;
        }
    }

    outputFile.close();

    return 0;
}

