/* kepler.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* Kepler Problem */

static int
rhs_kepler (double t, const double y[], double f[], void *p)
{
  (void) t, p;
  double r;

  r = sqrt (y[2] * y[2] + y[3] * y[3]);
  r = gsl_pow_3 (r);

  f[0] = -y[2] / r;
  f[1] = -y[3] / r;
  f[2] = y[0];
  f[3] = y[1];

  return GSL_SUCCESS;
}

static int
jac_kepler (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t, p;
  double r, r3, r5;

  r = sqrt (y[2] * y[2] + y[3] * y[3]);
  r3 = gsl_pow_3 (r);
  r5 = r3 * gsl_pow_2 (r);

  dfdy[0] = dfdy[1] = 0;
  dfdy[2] = -1.0 / r3 + 3.0 * y[2] * y[2] / r5;
  dfdy[3] = 3.0 * y[2] * y[3] / r5;
  dfdy[4] = dfdy[5] = 0;
  dfdy[6] = dfdy[3];
  dfdy[7] = -1.0 / r3 + 3.0 * y[3] * y[3] / r5;
  dfdy[8] = 1.0;
  dfdy[9] = dfdy[10] = dfdy[11] = dfdy[12] = 0;
  dfdy[13] = 1.0;
  dfdy[14] = dfdy[15] = 0;

  dfdt[0] = dfdt[1] = dfdt[2] = dfdt[3] = 0;

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_kepler = { rhs_kepler, jac_kepler, 4, 0 };

int
int_kepler (double t, const double y[], double I[], void *p)
{
  (void) t, p;

  I[0] = y[0] * y[0] + y[1] * y[1];
  I[0] /= 2;
  I[0] -= 1 / sqrt (y[2] * y[2] + y[3] * y[3]);

  I[1] = y[0] * y[3] - y[1] * y[2];

  return GSL_SUCCESS;
}

int
exact_kepler (const double t, double y[], const double y0[], void *p)
{
  (void) p;

  /* TODO: exact solution for arbitrary y0 & t */

  const double K = (y0[0] * y0[0] + y0[1] * y0[1]) / 2.0;
  const double U = -1.0 / sqrt (y0[2] * y0[2] + y0[3] * y0[3]);
  const double E = K + U;
  const double L = y0[1] * y0[2] - y0[0] * y0[3];
  const double e = sqrt (1 + 2 * E * L * L);
  const double a = -1 / 2.0 / E;
  const double xi = 4 * M_PI;
  /* const double t = sqrt (gsl_pow_3 (a)) * xi; */

  y[0] = sin (xi) / sqrt (a) / (e * cos (xi) - 1);
  y[1] = sqrt (1 - e * e) * cos (xi) / sqrt (a) / (1 - e * cos (xi));
  y[2] = a * (cos (xi) - e);
  y[3] = a * sqrt (1 - e * e) * sin (xi);

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_kepler = &sys_kepler;
