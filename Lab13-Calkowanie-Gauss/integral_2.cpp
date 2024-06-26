#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"
#include "/usr/include/gsl/gsl_eigen.h"
#include "/usr/include/gsl/gsl_errno.h"
#include "/usr/include/gsl/gsl_fft_complex.h"
#include "/usr/include/gsl/gsl_integration.h"
#include <math.h>

#define EPS 3.0e-11
#define MAXIT 10




double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
void nrerror(const char error_text[]) {
    // Implementation of nrerror function
    // This function can be used to throw an exception in C++.
    throw std::runtime_error(error_text);
}
void gaulag(std::vector<double>& x, std::vector<double>& w, int n, double alf) {
    int i, its, j;
    double ai;
    double p1, p2, p3, pp, z, z1;

    x.resize(n);
    w.resize(n);

    for (i = 0; i < n; i++) {
        if (i == 0) {
            z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
        } else if (i == 1) {
            z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
        } else {
            ai = i - 1;
            z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf /
                (1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * alf);
        }
        for (its = 0; its < MAXIT; its++) {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
            }
            pp = (n * p1 - (n + alf) * p2) / z;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS) break;
        }
        if (its >= MAXIT) nrerror("too many iterations in gaulag");
        x[i] = z;
        w[i] = -exp(gammln(alf + n) - gammln((double)n)) / (pp * n * p2);
    }
}


double integrand(double x, double k){
    return pow(x, k);
}

int main(void) {
    int n = 20;
    int k1 = 5;
    int k2 = 10;

    std::cout << "\n\n=====================[k=05]=====================\n\n";

    for (int i = 2; i <= n; ++i) {
        std::vector<double> x;
        std::vector<double> w;
        gaulag(x, w, i, 0);
        
        double integral = 0.0;
        for(int j=0;j<=i;j++){
            integral += w[j] * integrand(x[j], k1);
        }
        cout <<i<<" "<<integral<< endl;
    }

    std::cout << "\n\n=====================[k=10]=====================\n\n";

    for (int i = 2; i <= n; ++i) {
        std::vector<double> x;
        std::vector<double> w;
        gaulag(x, w, i, 0);
        
        double integral = 0.0;
        for(int j=0;j<=i;j++){
            integral += w[j] * integrand(x[j], k2);
        }
        cout <<i<<" "<<integral<< endl;
    }



    return 0;
}
