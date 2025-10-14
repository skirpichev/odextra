/* fermi.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_pow_int.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* Fermi–Pasta–Ulam problem */

static double alpha_fermi = 0.25;

static int
rhs_fermi (double t, const double y[], double f[], void *p)
{
  (void) t;

  const size_t n = 128, d = n / 2;

  double alpha = *(double*) p;

  f[0] = y[d + 1] - 2 * y[d] +
    alpha * (gsl_pow_2 (y[d + 1] - y[d]) - y[d] * y[d]);
  f[d] = y[0];

  f[d - 1] = y[n - 2] - 2 * y[n - 1] +
    alpha * (gsl_pow_2 (y[n - 1]) - gsl_pow_2 (y[n - 1] - y[n - 2]));
  f[n - 1] = y[d - 1];

  for (size_t i = 1; i < d - 1; i++)
    {
      f[i] = y[d + i + 1] + y[d + i - 1] - 2 * y[d + i] +
        alpha * (gsl_pow_2 (y[d + i + 1] - y[d + i]) -
                 gsl_pow_2 (y[d + i - 1] - y[d + i]));
      f[i + d] = y[i];
    }

  return GSL_SUCCESS;
}

static int
jac_fermi (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t;
  const size_t n = 128, d = n / 2;
  double alpha = *(double*) p;
  gsl_matrix_view dfdy_m = gsl_matrix_view_array (dfdy, n, n);
  gsl_matrix * m = &dfdy_m.matrix;

  gsl_matrix_set_zero (m);

  for (size_t i = 0; i < d; i++)
    {
      gsl_matrix_set (m, d + i, i, 1.0);

      dfdt[i] = dfdt[i + d] = 0;
    }

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_fermi = { rhs_fermi, jac_fermi, 128, &alpha_fermi };

int
int_fermi (double t, const double y[], double I[], void *p)
{
  (void) t;

  const size_t n = 128, d = n / 2;

  double alpha = *(double*) p;

  I[0] = y[0] * y[0] / 2 + y[d - 1] * y[d - 1] / 2;
  I[0] += gsl_pow_2 (y[n - 1] - y[d]) / 2;
  I[0] += alpha * gsl_pow_3 (y[n - 1] - y[d]) / 3;

  for (size_t i = 1; i < d - 1; i++)
    {
      I[0] += y[i] * y[i] / 2;

      I[0] += gsl_pow_2 (y[d + i] - y[d + i + 1]) / 2;
      I[0] += alpha * gsl_pow_3 (y[d + i] - y[d + i + 1]) / 3;
    }

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_fermi = &sys_fermi;
