#include <stdio.h>
#include <gsl/gsl_linalg.h>

int
main (void)
{
  double a_data[] = { -1, 2,  1, -1,
                      1, -3, -2, -1,
                      3, -1,  1, 4,
};

  double b_data[] = { 1.0, 2.0, 3.0 };

  gsl_matrix_view m
    = gsl_matrix_view_array (a_data, 3, 4);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, 3);

  gsl_vector *x = gsl_vector_alloc (4);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (4);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  printf ("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
  return 0;
}