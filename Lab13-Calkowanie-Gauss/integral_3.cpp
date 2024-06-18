#include <iostream>
#include <cmath>
#include <vector>

#define EPS 3.0e-11
#define MAXIT 10
#define PIM4 0.7511255444649425

void nrerror(const char error_text[]) {
    throw std::runtime_error(error_text);
}

void gauher(std::vector<double>& x, std::vector<double>& w, int n) {
    int i, its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n + 1) / 2;
    x.resize(n);
    w.resize(n);

    for (i = 0; i < m; i++) {
        if (i == 0) {
            z = sqrt(2.0 * n + 1) - 1.85575 * pow(2.0 * n + 1, -0.16667);
        } else if (i == 1) {
            z -= 1.14 * pow(n, 0.426) / z;
        } else if (i == 2) {
            z = 1.86 * z - 0.86 * x[0];
        } else if (i == 3) {
            z = 1.91 * z - 0.91 * x[1];
        } else {
            z = 2.0 * z - x[i - 2];
        }
        for (its = 0; its < MAXIT; its++) {
            p1 = PIM4;
            p2 = 0.0;
            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / j) * p2 - sqrt((j - 1.0) / j) * p3;
            }
            pp = sqrt(2.0 * n) * p2;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS) break;
        }
        if (its >= MAXIT) nrerror("too many iterations in gauher");
        x[i] = z;
        x[n - 1 - i] = -z;
        w[i] = 2.0 / (pp * pp);
        w[n - 1 - i] = w[i];
    }
}

double integrand(double x, double y) {
    return sin(x) * sin(x) * sin(y) * sin(y) * sin(y) * sin(y);
}

int main() {
    int max_n = 15;
    double exact = 0.1919832644;

    std::cout << "n\tApprox Integral\tExact Integral\tError" << std::endl;

    for (int n = 2; n <= max_n; ++n) {
        std::vector<double> x, w;
        gauher(x, w, n);

        double integral = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                integral += w[i] * w[j] * integrand(x[i], x[j]);
            }
        }

        double error = fabs(integral - exact);
        std::cout << n << "\t" << integral << "\t" << exact << "\t" << error << std::endl;
    }


    return 0;
}
