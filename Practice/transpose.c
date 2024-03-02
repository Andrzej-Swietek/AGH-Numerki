#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

int main() {
    gsl_matrix *A = gsl_matrix_alloc(2, 3); // A: 2x3 matrix
    gsl_matrix *AT = gsl_matrix_alloc(3, 2); // AT: Transpose of A
    gsl_vector *v = gsl_vector_alloc(3); // v: Vector of size 3

    // Initialize matrix A and vector v with some values
    gsl_matrix_set_all(A, 1.0);
    gsl_vector_set_all(v, 2.0);

    // Transpose matrix A: AT = A^T
    gsl_matrix_transpose_memcpy(AT, A);

    // Print the transpose matrix AT
    printf("Transpose of matrix A:\n");
    gsl_matrix_fprintf(stdout, AT, "%g");

    // Perform vector-matrix multiplication: result = A * v
    gsl_vector *result = gsl_vector_alloc(2); // Resultant vector
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, v, 0.0, result);

    // Print the resultant vector
    printf("\nResult of vector-matrix multiplication:\n");
    gsl_vector_fprintf(stdout, result, "%g");

    // Free allocated memory
    gsl_matrix_free(A);
    gsl_matrix_free(AT);
    gsl_vector_free(v);
    gsl_vector_free(result);

    return 0;
}
