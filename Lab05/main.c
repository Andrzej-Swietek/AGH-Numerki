#include <stdio.h>
#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"
#include "/usr/include/gsl/gsl_eigen.h"


double rho(double x, double alpha) {
    return 1 + 4*alpha*x*x;
}

int delta_kroneckera(int i, int j){
    return i == j ? 1 : 0;
}


int main() {
    const int n = 200;
    const int L = 10; // dl struny
    const int N = 1;

    gsl_vector *eval = gsl_vector_calloc(n);
    gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);
    
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_matrix *B = gsl_matrix_calloc(n, n);
    gsl_matrix *evec = gsl_matrix_calloc(n, n); // eigen vectors

    double deltaX = 10.0/(200.0 - 1.0); // L/(n-1)
    FILE *fp_alpha = fopen("data/alpha_values.txt", "w");
    FILE *fp_eigenvalues = fopen("data/eigenvalues.txt", "w");

    fprintf(fp_alpha, "alpha\n");
    fprintf(fp_eigenvalues, "alpha\teigenvalue_1\teigenvalue_2\teigenvalue_3\teigenvalue_4\teigenvalue_5\teigenvalue_6\n");
    
    for(int alpha = 0; alpha < 100; alpha+=2) {  // Pojedyncy problem


        // Wypełnienie Macierzy A i B
        for(int i =0; i < n; i++){
            for(int j = 0; j < n; j++) {
                double value = (-delta_kroneckera(i, j+1) + 2*delta_kroneckera(i,j) - delta_kroneckera(i,j-1))/(deltaX*deltaX);
                gsl_matrix_set(A,i,j,value);
                double beta_i_j = rho(i, alpha)/N * delta_kroneckera(i,j);
                gsl_matrix_set(B,i,j, beta_i_j);
                // printf("%lf ", value);
            }
            // printf("\n");
        }

        gsl_eigen_gensymmv(A, B, eval, evec, w);
        gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        // Zapisanie wartości pierwiastków do pliku
        fprintf(fp_alpha, "%d\n", alpha);
        for (int i = 0; i < 6; i++) {
            double eigenvalue = sqrt(gsl_vector_get(eval, i));
            fprintf(fp_eigenvalues, "%d\t%f", alpha, eigenvalue);
        }
        fprintf(fp_eigenvalues, "\n");
    
        if (alpha == 0 || alpha == 98) {
            char filename[50];
            sprintf(filename, "data/eigenvectors_alpha_%d.txt", alpha);
            FILE *fp_eigenvectors = fopen(filename, "w");
            for (int i = 0; i < 6; i++) {
                fprintf(fp_eigenvectors, "Eigenvalue %d\n", i+1);
                for (int j = 0; j < n; j++) {
                    double eigenvector_element = gsl_matrix_get(evec, j, i);
                    fprintf(fp_eigenvectors, "%f\n", eigenvector_element);
                }
            }
            fclose(fp_eigenvectors);
        }
    }

    fclose(fp_alpha);
    fclose(fp_eigenvalues);

    gsl_vector_free(eval);
    gsl_eigen_gensymmv_free(w);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);
    return 0;
}