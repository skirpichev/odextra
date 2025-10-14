/* linear.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_math.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* Harmonic Oscillator */

static int
rhs_linear (double t, const double y[], double f[], void *p)
{
  (void) t, p;

  f[0] = -y[2];
  f[1] = -y[3];
  f[2] = y[0];
  f[3] = y[1];

  return GSL_SUCCESS;
}

static int
jac_linear (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t, p;

  dfdy[0] = dfdy[1] = 0;
  dfdy[2] = -1.0;
  dfdy[3] = dfdy[4] = dfdy[5] = dfdy[6] = 0;
  dfdy[7] = -1.0;
  dfdy[8] = 1.0;
  dfdy[9] = dfdy[10] = dfdy[11] = dfdy[12] = 0;
  dfdy[13] = 1.0;
  dfdy[14] = dfdy[15] = 0;

  dfdt[0] = dfdt[1] = dfdt[2] = dfdt[3] = 0;

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_linear = { rhs_linear, jac_linear, 4, 0 };

int
int_linear (const double t, const double y[], double I[], void *p)
{
  (void) t, p;

  I[0] = 0;
  for (size_t i = 0; i < 4; i++)
    {
      I[0] = y[i] * y[i];
    }
  I[0] /= 2;

  return GSL_SUCCESS;
}

int
exact_linear (const double t, double y[], const double y0[], void *p)
{
  (void) p;

  double c = cos (t), s = sin (t);

  y[0] = y0[0] * c - y0[2] * s;
  y[1] = y0[1] * c - y0[3] * s;
  y[2] = y0[0] * s + y0[2] * c;
  y[3] = y0[1] * s + y0[3] * c;

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_linear = &sys_linear;
