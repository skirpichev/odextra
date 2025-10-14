/* vortex.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_pow_int.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* The n-vortex problem */

static double gamma_vortex[4] = {1, 2, 3, 4};

static int
rhs_vortex (double t, const double y[], double f[], void *p)
{
  (void) t;

  const size_t n = 8, d = n / 2;

  double* gamma = (double*) p;

  for (size_t i = 0; i < d; i++)
    {
      double zi = y[i] / gamma[i];

      f[i] = 0;
      f[i + d] = 0;

      for (size_t j = 0; j < d; j++)
        {
          if (i != j)
            {
              double zj = zi - y[j] / gamma[j];
              double r = zj * zj + gsl_pow_2 (y[i + d] - y[j + d]);

              f[i] -= gamma[i] * gamma[j] * (y[i + d] - y[j + d]) / r / 2.0;
              f[i + d] += gamma[j] * zj / r / 2.0;
            }
        }
    }

  return GSL_SUCCESS;
}

static int
jac_vortex (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t;

  const size_t n = 8, d = n / 2;
  double* gamma = (double*) p;
  gsl_matrix_view dfdy_m = gsl_matrix_view_array (dfdy, n, n);
  gsl_matrix * m = &dfdy_m.matrix;

  gsl_matrix_set_zero (m);

  for (size_t i = 0; i < d; i++)
    {
      double zi = y[i] / gamma[i];
      double diag[4] = {0, 0, 0, 0};

      for (size_t j = 0; j < d; j++)
        {
          double zj = zi - y[j] / gamma[j];
          double r2 = gsl_pow_2 (zj * zj + gsl_pow_2 (y[i + d] - y[j + d]));

          if (i != j)
            {
              double v;

              v = gamma[i] * (y[i + d] - y[j + d]) * zj / r2;
              gsl_matrix_set (m, i, j, v);
              diag[0] -= v;

              v = gamma[i] * gamma[j] * (0.5 - gsl_pow_2 (y[i + d] - y[j + d]) / r2) / r2;
              gsl_matrix_set (m, i, j + d, v);
              diag[1] -= v;

              v = (0.5 - zj * zj / r2) / r2;
              gsl_matrix_set (m, i + d, j, v);
              diag[2] -= v;

              v = -gamma[j] * zj * (y[i + d] - y[j + d]) / gsl_pow_2 (r2);
              gsl_matrix_set (m, i + d, j + d, v);
              diag[3] -= v;
            }
        }

      gsl_matrix_set (m, i, i, diag[0]);
      gsl_matrix_set (m, i, i + d, diag[1]);
      gsl_matrix_set (m, i + d, i, diag[2]);
      gsl_matrix_set (m, i + d, i + d, diag[3]);

      dfdt[i] = dfdt[i + d] = 0;
    }

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_vortex = { rhs_vortex, jac_vortex, 8, gamma_vortex };

int
int_vortex (double t, const double y[], double I[], void *p)
{
  (void) t;

  const size_t n = 8, d = n / 2;

  double* gamma = (double*) p;

  I[0] = I[1] = I[2] = I[3] = 0;

  for (size_t i = 0; i < d; i++)
    {
      double zi = y[i] / gamma[i];

      I[1] += gamma[i] * (y[i + d] * y[i + d] + zi * zi); /* L^2 */
      I[2] += y[i]; /* P */
      I[3] += gamma[i] * y[i + d]; /* Q */

      for (size_t j = 0; j < d; j++)
        {
          if (i != j)
            {
              double zj = zi - y[j] / gamma[j];
              double r = zj * zj + gsl_pow_2 (y[i + d] - y[j + d]);

              I[0] -= gamma[i] * gamma[j] * log (r) / 4.0; /* H */
            }
        }
    }

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_vortex = &sys_vortex;
