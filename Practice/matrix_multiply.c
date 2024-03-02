#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int main() {
    gsl_matrix *A = gsl_matrix_alloc(2, 3); // A: 2x3 matrix
    gsl_matrix *B = gsl_matrix_alloc(3, 2); // B: 3x2 matrix
    gsl_matrix *C = gsl_matrix_alloc(2, 2); // C: Resultant 2x2 matrix

    // Initialize matrices A and B with some values
    gsl_matrix_set_all(A, 1.0);
    gsl_matrix_set_all(B, 2.0);

    // Perform matrix multiplication: C = A * B
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);

    // Print the result matrix C
    printf("Resultant matrix C:\n");
    gsl_matrix_fprintf(stdout, C, "%g");

    // Free allocated memory
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(C);

    return 0;
}
