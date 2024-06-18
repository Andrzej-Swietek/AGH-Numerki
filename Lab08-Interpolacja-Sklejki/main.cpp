#include<iostream>
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

double f1(double x) {
    return 1.0/(1.0 + x*x);
}

double f2(double x) {
    return cos(2*x);
}

double phi(double x, double xi_minus_2, double xi_minus_1, double xi, double xi_plus_1, double xi_plus_2, double h) {
    if( x>= xi_minus_2 && x < xi_minus_1 ) return 1/pow(h,3) * pow(x - xi_minus_2, 3);
    else if( x >= xi_minus_1 && x < xi )
        return 1/pow(h,3) * ( pow(h,3) + 3*pow(h,2)*(x - xi_minus_1) + 3*h* pow((x - xi_minus_1),2) - 3* pow((x - xi_minus_1),3) );
    else if( x >= xi && x < xi_plus_1)
        return 1/pow(h,3) * ( pow(h,3) + 3*pow(h,2)*(xi_plus_1 - x) + 3*h* pow((xi_plus_1 - x),2) - 3* pow((xi_plus_1 - x),3) );
    else if ( x >= xi_plus_1 && x < xi_plus_2)
        return 1/pow(h,3) * pow((xi_plus_2 - x) ,3);
    return  0.0;
}


double interpolate(double x, const std::vector<double>& c, const std::vector<double>& xx, double h) {
    double sum = 0.0;
    int n = c.size(); // N+1 - po spuszowaniu C_0 i C_n+1
    cout << endl << endl;
    for (int i = 0; i < n ; ++i) {
        sum += c[i] * phi(x, xx[i], xx[i+1], xx[i+2], xx[i+3], xx[i+4], h);
    }
    return sum;
}


std::vector<double> convert_gsl_to_vector(gsl_vector* vector_gsl, int N){
        std::vector<double> vect(N);
        for (int i = 0; i < N; ++i) {
            vect[i] = gsl_vector_get(vector_gsl, i);
        }
        return vect;
}


double dfdx(std::function<double(double)>f, double x) {
    const double Delta_X = 0.01;
    return ( f(x+Delta_X) - f(x-Delta_X) ) / (2*Delta_X);
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

void saveNodesToFile(const vector<double>& x, const vector<double>& y, int num_nodes) {
    ofstream file("data/nodes_xi_" + to_string(num_nodes) + ".txt");
    file << "x" << " " << "y\n";
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << " " << y[i] << endl;
    }
    file.close();
}



int main() {
    const int num_nodes[] = {5,6, 7,10, 20};
    const double x_min = -5.0; const double x_max = 5.0;

    double alpha = dfdx(f1, x_min);
    double beta  = dfdx(f1, x_max);

    for( int N : num_nodes) {
        cout << "===============[ " << N << " ]===============\n";
        double h = (x_max - x_min) / (N - 1);
        std::vector<double> xx(N + 6); // z dołożonymi wezlami
        std::vector<double> xw(N + 4);
        std::vector<double> yw(N + 4);

        cout << "H:" << h <<"\n";

       // WEZLY
        for (int i = -2; i <= (N + 3); ++i)
            xx[i+2] = x_min + h * (i-1);


        for (int i = 1; i <= (N + 4); ++i)
            xw[i-1] = xx[i+1];

        cout << "---------------------------------------- \n";
        for (int i = 0; i <= (N + 3); ++i){
            yw[i] = f1(xw[i]);
            cout << xw[i] << " " << yw[i] << " " << f1(xw[i]) << endl;
        }
        cout << "---------------------------------------- \n";

        saveNodesToFile(xw, yw, N);
        

        // ROZWIAZNANIE UKL ROWNAN
        gsl_matrix *A = gsl_matrix_alloc(N, N);
        gsl_vector *b = gsl_vector_alloc(N);
        gsl_vector *c_gsl = gsl_vector_alloc(N);

        // Wypełnienie macierzy A i wektora b
        for (int i = 0; i < N; ++i) {
            if (i == 0 || i == N-1) {
                gsl_matrix_set(A, i, i, 4.0);
                if (i == 0) {
                    gsl_matrix_set(A, i, i + 1, 2.0);
                } else {
                    gsl_matrix_set(A, i, i - 1, 2.0);
                }
                gsl_vector_set(
                        b, i, i == 0 ? f1(xw[i+1]) + alpha*h / 3 : f1(xw[i+1]) - beta*h / 3
                );
            } else {
                gsl_matrix_set(A, i, i, 4.0);
                gsl_matrix_set(A, i, i - 1, 1.0);
                gsl_matrix_set(A, i, i + 1, 1.0);
                gsl_vector_set(b, i, f1(xw[i+1]));
            }
        }
        
//        print_gsl_matrix(A);

        std::cout << "\n[Vector x]: ";
        for ( double value : xx ) cout << value << " ";
        std::cout << "\n";

        std::cout << "\n[Vector xw]: ";
        for ( double value : xw ) cout << value << " ";
        std::cout << "\n";

        std::cout << "\n[Vector B]: ";
        for (int i = 0; i < N; i++) {
            std::cout << gsl_vector_get(b, i) << " ";
        }
        std::cout << "\n";

        // Rozwiązanie układu równań
        gsl_linalg_HH_solve(A, b, c_gsl);
//        gsl_permutation *p = gsl_permutation_alloc(N);
//        int signum;
//        gsl_linalg_LU_decomp(A, p, &signum);
//        gsl_linalg_LU_solve(A, p, b, c_gsl);


        std::cout << "\n[Vector C GSL]: ";
        for (int i = 0; i < N; i++) {
            std::cout << gsl_vector_get(c_gsl, i) << " ";
        }
        std::cout << "\n";

        vector<double> c = convert_gsl_to_vector(c_gsl, N);
        std::cout << "[Vector C]: ";
        for (double value : c)
            std::cout << value << " ";
        std::cout << "\n";

          
        double step = 0.1;
        ofstream file("data/interpolation_results_" + to_string(N) + ".txt");
        file << "x y\n";

        double c0 = c[1] - h/3*alpha;
        double cn_plus_1 = c[ c.size() -2] + h/3*beta;

        c.insert(c.begin(), c0);
        c.push_back(cn_plus_1);

        std::cout << "[Vector C_new]: ";
        for (double value : c)
            std::cout << value << " ";
        std::cout << "\n";

        std::cout << "Number of nodes (n): " << N << std::endl;
        for (double x = x_min; x <= x_max; x += step) {
            double s_x = interpolate(x, c, xx, h);
//            cout << x << " " << s_x << "\n";
            file << x << " " << s_x << "\n";
        }
        file.close();
        std::cout << std::endl;

 
        gsl_matrix_free(A);
        gsl_vector_free(b);
        gsl_vector_free(c_gsl);
    }

    return 0;
}
