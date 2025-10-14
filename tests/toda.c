/* toda.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* Toda Lattice */

static int
rhs_toda (double t, const double y[], double f[], void *p)
{
  (void) t, p;

  const size_t n = 20, d = n / 2;

  for (size_t i = 0; i < d; i++)
    {
      double r,l;

      if (i < d - 1)
        r = y[d + i + 1];
      else
        r = y[d];

      if (i > 0)
        l = y[d + i - 1];
      else
        l = y[n - 1];

      f[i] = exp (l - y[d + i]) - exp (y[d + i] - r);
      f[i + d] = y[i];
    }

  return GSL_SUCCESS;
}

static int
jac_toda (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t, p;
  const size_t n = 20, d = n / 2;
  gsl_matrix_view dfdy_m = gsl_matrix_view_array (dfdy, n, n);
  gsl_matrix * m = &dfdy_m.matrix;

  gsl_matrix_set_zero (m);

  for (size_t i = 0; i < d; i++)
    {
      double l, r, v;

      if (i < d - 1)
        r = y[d + i + 1];
      else
        r = y[d];

      if (i > 0)
        l = y[d + i - 1];
      else
        l = y[n - 1];

      v = -exp (y[d + i] - r) - exp (l - y[d + i]);
      gsl_matrix_set (m, i + d, i, v);
      gsl_matrix_set (m, i, i + d, v);

      v = exp (y[d + i] - r);
      gsl_matrix_set (m, i + d, i + 1, v);
      gsl_matrix_set (m, i, i + d + 1, v);

      v = exp (l - y[d + i]);
      gsl_matrix_set (m, i + d, i - 1, v);
      gsl_matrix_set (m, i, i + d - 1, v);

      gsl_matrix_set (m, d + i, d + i, 1.0);

      dfdt[i] = dfdt[i + d] = 0;
    }

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_toda = { rhs_toda, jac_toda, 20, 0 };

int
int_toda (double t, const double y[], double I[], void *p)
{
  (void) t, p;

  const size_t n = 20, d = n / 2;

  /* Find first integrals of motion (Lax pair) */

  gsl_matrix * L = gsl_matrix_alloc (d, d);

  gsl_matrix_set_zero (L);

  for (size_t i = 0; i < d; i++)
    {
      gsl_matrix_set (L, i, i, -y[i] / 2);

      if (i < d - 1)
        {
          double v = exp ((y[d + i] - y[d + i + 1]) / 2) / 2;
          gsl_matrix_set (L, i, i + 1, v);
          gsl_matrix_set (L, i + 1, i, v);
        }
    }

  {
    double v = exp ((y[n - 1] - y[d]) / 2) / 2;
    gsl_matrix_set (L, 0, d - 1, v);
    gsl_matrix_set (L, d - 1, 0, v);
  }

  gsl_vector i = gsl_vector_view_array (I, d).vector;
  gsl_matrix * tmp = gsl_matrix_alloc (d, d);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (d);
  gsl_eigen_symmv (L, &i, tmp, w);
  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (&i, tmp, GSL_EIGEN_SORT_VAL_ASC);

  gsl_matrix_free (tmp);
  gsl_matrix_free (L);

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_toda = &sys_toda;
