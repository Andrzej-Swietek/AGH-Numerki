#include <iostream>
#include <cmath>
#include <vector>

void gauleg(double x1, double x2, std::vector<double>& x, std::vector<double>& w, int n) {
    double eps = 1.0e-14;
    int m = (n + 1) / 2;
    double xm = 0.5 * (x2 + x1);
    double xl = 0.5 * (x2 - x1);

    x.resize(n);
    w.resize(n);

    for (int i = 0; i < m; i++) {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1;
	double pp;
        do {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 0; j < n; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > eps);

        x[i] = xm - xl * z;
        x[n - 1 - i] = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

double integrand(double x) {
    return x / (4 * x * x + 1);
}

int main() {
    double a = 0;
    double b = 2;
    int max_n = 20;

    std::cout << "n\tApprox Integral\tExact Integral\tError" << std::endl;
    //double exact = 1/8 * log(16.0 / 3.0); // derived from the provided analytical solution
    double exact = 0.3542;
    for (int n = 2; n <= max_n; ++n) {
        std::vector<double> x, w;
        gauleg(a, b, x, w, n);

        double integral = 0.0;
        for (int i = 0; i < n; i++) {
            integral += w[i] * integrand(x[i]);
        }

        double error = fabs(integral - exact);
        std::cout << n << "\t" << integral << "\t" << exact << "\t" << error << std::endl;
    }

    return 0;
}