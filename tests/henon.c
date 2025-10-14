/* henon.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

#include "odeiv_util.h"
#include "odextra_test.h"

/* Henon-Heiles Hamiltonian */

static int
rhs_henon (double t, const double y[], double f[], void *p)
{
  (void) t, p;

  f[0] = -y[2] - 2.0 * y[2] * y[3];
  f[1] = -y[3] - y[2] * y[2] + y[3] * y[3];
  f[2] = y[0];
  f[3] = y[1];

  return GSL_SUCCESS;
}

static int
jac_henon (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t, p;

  dfdy[0] = dfdy[1] = 0;
  dfdy[2] = -1.0 - 2.0 * y[3];
  dfdy[3] = -2.0 * y[2];
  dfdy[4] = dfdy[5] = 0;
  dfdy[6] = dfdy[3];
  dfdy[7] = -1.0 + 2.0 * y[3];
  dfdy[8] = 1.0;
  dfdy[9] = dfdy[10] = dfdy[11] = dfdy[12] = 0;
  dfdy[13] = 1.0;
  dfdy[14] = dfdy[15] = 0;

  dfdt[0] = dfdt[1] = dfdt[2] = dfdt[3] = 0;

  return GSL_SUCCESS;
}

static const gsl_odeiv_system sys_henon = { rhs_henon, jac_henon, 4, 0 };

int
int_henon (double t, const double y[], double I[], void *p)
{
  (void) t, p;

  I[0] = 0.0;
  for (size_t i = 0; i < 4; i++)
    {
      I[0] += y[i] * y[i];
    }
  I[0] /= 2;

  I[0] += (y[2] * y[2] - y[3] * y[3] / 3.0) * y[3];

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_henon = &sys_henon;
