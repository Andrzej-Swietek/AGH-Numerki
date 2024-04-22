//
// Created by andrzej on 22.04.24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

// Funkcja obliczająca węzły Czebyszewa
vector<double> chebyshevNodes(int num_nodes, double xa, double xb) {
    vector<double> nodes(num_nodes + 1);
    for (int i = 0; i <= num_nodes; ++i) {
        nodes[i] = (xa + xb) / 2.0 + 0.5 * (xb - xa) * cos((2.0 * i + 1.0) * M_PI / (2.0 * num_nodes + 2.0));
    }
    return nodes;
}

// Implementacja interpolacji Lagrange'a
double lagrangeInterpolation(const vector<double>& x, const vector<double>& y, double xi) {
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double term = y[i];
        for (size_t j = 0; j < x.size(); ++j) {
            if (j != i) {
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

// Funkcja zapisująca węzły do pliku
void saveNodesToFile(const vector<double>& x, const vector<double>& y, int num_nodes) {
    ofstream file("data/nodes_xi_" + to_string(num_nodes) + "_czebyszew.txt");
    file << "x" << " " << "y\n";
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << " " << y[i] << endl;
    }
    file.close();
}

// Funkcja zapisująca wyniki interpolacji do pliku
void saveInterpolationResultsToFile(const vector<double>& x_values, const vector<double>& y_values, const vector<double>& interpolated_values, int num_nodes) {
    ofstream file("data/interpolation_results_" + to_string(num_nodes) + "_czebyszew.txt");
    file << "x y interpolated\n";
    for (size_t i = 0; i < x_values.size(); ++i) {
        file << x_values[i] << " " << y_values[i] << " " << interpolated_values[i] << endl;
    }
    file.close();
}

int main() {
    const int num_nodes[] = {5, 10, 15};
    const double xa = -3.0;
    const double xb = 8.0;

    // Funkcja y(x) = x / (1 + x^2)
    auto y_func = [](double x) { return x / (1 + x * x); };

    // Wykresy
    const int num_points = 200;
    vector<double> x_values(num_points), y_values(num_points);

    double dx = (xb - xa) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        x_values[i] = xa + i * dx;
        y_values[i] = y_func(x_values[i]);
    }

    for (int n : num_nodes) {
        // Węzły równomiernie rozłożone
        vector<double> x_uniform(n + 1), y_uniform(n + 1);
        for (int i = 0; i <= n; ++i) {
            x_uniform[i] = xa + i * (xb - xa) / n;
            y_uniform[i] = y_func(x_uniform[i]);
        }
        saveNodesToFile(x_uniform, y_uniform, n);

        // Węzły Czebyszewa
        vector<double> x_chebyshev = chebyshevNodes(n, xa, xb);
        vector<double> y_chebyshev(n + 1);
        for (int i = 0; i <= n; ++i) {
            y_chebyshev[i] = y_func(x_chebyshev[i]);
        }
        saveNodesToFile(x_chebyshev, y_chebyshev, n);

        // Interpolacja Lagrange'a dla węzłów Czebyszewa
        vector<double> interpolated_chebyshev_values(num_points);
        for (int i = 0; i < num_points; ++i) {
            interpolated_chebyshev_values[i] = lagrangeInterpolation(x_chebyshev, y_chebyshev, x_values[i]);
        }
        saveInterpolationResultsToFile(x_values, y_values, interpolated_chebyshev_values, n);
    }

    return 0;
}
