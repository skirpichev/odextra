/* pendulum.c */

#include <config.h>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_sf_ellint.h>

#include "odextra_test.h"
#include "odeiv_util.h"

/* Pendulum */

static int
rhs_pendulum (double t, const double y[], double f[], void *p)
{
  (void) t, p;

  f[0] = -sin (y[1]);
  f[1] = y[0];

  return GSL_SUCCESS;
}

static int
jac_pendulum (double t, const double y[], double *dfdy, double dfdt[], void *p)
{
  (void) t, p;

  dfdy[0] = 0;
  dfdy[1] = -cos (y[1]);
  dfdy[2] = 1.0;
  dfdy[3] = 0;

  dfdt[0] = dfdt[1] = 0;

  return GSL_SUCCESS;
}

static gsl_odeiv_system sys_pendulum = { rhs_pendulum, jac_pendulum, 2, 0 };

int
int_pendulum (double t, const double y[], double I[], void *p)
{
  (void) t, p;

  I[0] = y[0] * y[0] / 2 - cos(y[1]);

  return GSL_SUCCESS;
}

int
exact_pendulum (const double t, double y[], const double y0[], void *p)
{
  (void) p;

  double k, sn, cn, dn, H;
  gsl_sf_result t0;

  int_pendulum (0, y0, &H, 0);

  k = sqrt ((1 + H) / 2);

  if (k < 1)
    {
      int s;

      H = -asin (sin (y0[1] / 2) / k);
      s = gsl_sf_ellint_F_e (H, k, GSL_PREC_DOUBLE, &t0);

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      s = gsl_sf_elljac_e (t - t0.val, k * k, &sn, &cn, &dn);

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      y[0] = 2.0 * k * cn;
      y[1] = 2.0 * asin (k * sn);
    }
  else if (k > 1)
    {
      return GSL_FAILURE; /* TODO */
    }
  else /* k == 1 */
    {
      return GSL_FAILURE; /* TODO */
    }

  return GSL_SUCCESS;
}

const gsl_odeiv_system *rhs_func_pendulum = &sys_pendulum;
