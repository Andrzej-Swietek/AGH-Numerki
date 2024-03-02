#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int main() {
    gsl_matrix *A = gsl_matrix_alloc(3, 3);
    gsl_matrix_set(A, 0, 0, 1.0);   gsl_matrix_set(A, 0, 1, 2.0);   gsl_matrix_set(A, 0, 2, 3.0);
    gsl_matrix_set(A, 1, 0, 4.0);   gsl_matrix_set(A, 1, 1, 5.0);   gsl_matrix_set(A, 1, 2, 6.0);
    gsl_matrix_set(A, 2, 0, 7.0);   gsl_matrix_set(A, 2, 1, 8.0);   gsl_matrix_set(A, 2, 2, 10.0);

    gsl_matrix *inverse = gsl_matrix_alloc(3, 3);

    gsl_permutation *p = gsl_permutation_alloc(3);
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);        // LU decomposition
    gsl_linalg_LU_invert(A, p, inverse);    // Compute the inverse

    printf("Inverse of A:\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%0.6f\t", gsl_matrix_get(inverse, i, j));
        }
        printf("\n");
    }

    gsl_matrix_free(A);
    gsl_matrix_free(inverse);
    gsl_permutation_free(p);

    return 0;
}
