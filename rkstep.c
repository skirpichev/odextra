/* rkstep.c */

/* Generic Runge-Kutta integrator. */

/* Error estimation by step doubling, see eg. Ascher, U.M., Petzold,
   L.R., Computer methods for ordinary differential and
   differential-algebraic equations, SIAM, Philadelphia, 1998. */

#include <config.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

rk_state_t*
rk_init (size_t dim, const double a[], const double b[],
         const double c[], const size_t s)
{
  rk_state_t *state = (rk_state_t *) malloc (sizeof (rk_state_t));

  if (state == NULL)
    {
      goto cleanup;
    }

  state->y0 = (double *) malloc (dim * sizeof (double));
  state->y_onestep = (double *) malloc (dim * sizeof (double));
  state->y0_orig = (double *) malloc (dim * sizeof (double));
  state->fY = (double *) malloc (s * dim * sizeof (double));
  state->fY_tmp = (double *) malloc (s * dim * sizeof (double));

  if (!state->y0 || !state->y_onestep || !state->y0_orig ||
      !state->fY || !state->fY_tmp)
    {
      goto cleanup;
    }

  state->s = s;

  state->a = (double *) malloc (s * s * sizeof (double));
  state->b = (double *) malloc (s * sizeof (double));
  state->c = (double *) malloc (s * sizeof (double));

  if (!state->a || !state->b || !state->c)
    {
      goto cleanup;
    }

  cblas_dcopy (s * s, a, 1, state->a, 1);
  cblas_dcopy (s, b, 1, state->b, 1);
  cblas_dcopy (s, c, 1, state->c, 1);

  return state;

 cleanup:
  free (state->fY);
  free (state->c);
  free (state->b);
  free (state->a);
  free (state->y0_orig);
  free (state->y_onestep);
  free (state->y0);
  free (state);

  GSL_ERROR_NULL ("failed to allocate space for rk_state_t",  GSL_ENOMEM);
}

static int
rkstep (double *y, const double y0[], double Z[], double Z_tmp[],
        const double h, const double t, const gsl_odeiv_system *sys,
        const double a[], const double b[], const double c[], const size_t s)
{
  int ret = GSL_FAILURE;
  size_t dim = sys->dimension;

  /* Solve equations: Z_i = h f(y0 + sum j=1..s a_{ij} Z_j) */

  const size_t maxiter = 20;
  size_t iter = 1;
  const double tol = GSL_DBL_EPSILON;
  double norm;

  do
    {
      cblas_dcopy (dim * s, Z, 1, Z_tmp, 1);

      norm = 0; /* reset norm */

      for (size_t j = 0; j < s; j++)
        {
          cblas_dcopy (dim, y0, 1, Z + j * dim, 1);
          cblas_dgemv (CblasRowMajor, CblasTrans, s, dim, 1, Z_tmp, dim, a + j * s, 1, 1, Z + j * dim, 1);

          ret = GSL_ODEIV_FN_EVAL (sys, t + c[j] * h, Z + j * dim, y);

          if (ret != GSL_SUCCESS)
            {
              return ret;
            }

          cblas_dscal (dim, h, y, 1);
          cblas_dcopy (dim, y, 1, Z + j * dim, 1);

          cblas_daxpy (dim, -1, Z_tmp + j * dim, 1, y, 1);
          norm = GSL_MAX (norm, cblas_dasum (dim, y, 1));
        }

      iter++;

      if (iter > maxiter)
        {
          GSL_ERROR ("reached iteration limit for "
                     "fixed-point equation solver", GSL_EMAXITER);
        }
    }
  while (norm > tol);

  /* y' = y + h sum i=1..s b_i f(y + Z_i) */

  cblas_dcopy (dim, y0, 1, y, 1);
  cblas_dgemv (CblasRowMajor, CblasTrans, s, dim, 1, Z, dim, b, 1, 1, y, 1);

  return ret;
}

int
rk_apply (void *vstate, size_t dim, double t, double h, double y[],
          double yerr[], const double dydt_in[],
          double dydt_out[], const gsl_odeiv_system * sys)
{
  int ret = GSL_FAILURE;
  rk_state_t *state = (rk_state_t *) vstate;

  double *y0 = state->y0;
  double *y_onestep = state->y_onestep;
  double *y0_orig = state->y0_orig;

  /* Error estimation is done by step doubling procedure */

  /* initialization step */

  cblas_dcopy (dim, y, 1, y0, 1);
  cblas_dcopy (dim, y, 1, y0_orig, 1); /* Save initial values for possible failures */

  if (dydt_in != NULL)
    {
      cblas_dcopy (dim, dydt_in, 1, state->fY, 1);
    }
  else
    {
      ret = GSL_ODEIV_FN_EVAL (sys, t, y, state->fY);

      if (ret != GSL_SUCCESS)
        {
          goto cleanup;
        }
    }

  for (size_t j = 1; j < state->s; j++)
    {
      cblas_dcopy (dim, state->fY, 1, state->fY + dim * j, 1);
    }

  for (size_t j = 0; j < state->s; j++)
    {
      cblas_dscal (dim, h * state->c[j], state->fY + j * dim, 1);
    }

  /* First traverse h with one step (save to y_onestep) */

  cblas_dcopy (dim, y, 1, y_onestep, 1);

  ret = rkstep (y_onestep, state->y0,
                state->fY, state->fY_tmp, h, t, sys,
                state->a, state->b, state->c, state->s);

  if (ret != GSL_SUCCESS)
    {
      goto cleanup;
    }

  /* Then with two steps with half step length (save to y) */

  ret = rkstep (y, state->y0,
                state->fY, state->fY_tmp, h / 2, t, sys,
                state->a, state->b, state->c, state->s);

  if (ret != GSL_SUCCESS)
    {
      goto cleanup;
    }

  cblas_dcopy (dim, y, 1, y0, 1);

  ret = rkstep (y, state->y0,
                state->fY, state->fY_tmp, h / 2, t + h / 2, sys,
                state->a, state->b, state->c, state->s);

  if (ret != GSL_SUCCESS)
    {
      goto cleanup;
    }

  /* Derivatives at output */

  if (dydt_out != NULL)
    {
      ret = GSL_ODEIV_FN_EVAL (sys, t + h, y, dydt_out);

      if (ret != GSL_SUCCESS)
        {
          goto cleanup;
        }
    }

  /* Error estimation */

  /* Denominator in step doubling error equation
   * yerr = 8 * 0.5 * | y(onestep) - y(twosteps) | / (2^order - 1)
   */

  cblas_dcopy (dim, y_onestep, 1, yerr, 1);
  cblas_daxpy (dim, -1, y, 1, yerr, 1);
  cblas_dscal (dim, 4.0 / ((1 << state->s) - 1), yerr, 1);

  return GSL_SUCCESS;

cleanup:
  cblas_dcopy (dim, y0_orig, 1, y, 1); /* Restore original y vector */

  return ret;
}

void
rk_free (void *vstate)
{
  rk_state_t *state = (rk_state_t *) vstate;

  free (state->y_onestep);
  free (state->y0_orig);
  free (state->y0);
  free (state->a);
  free (state->b);
  free (state->c);
  free (state->fY);
  free (state->fY_tmp);
  free (state);
}

unsigned int
rk_order (void *vstate)
{
  rk_state_t *state = (rk_state_t *) vstate;

  return state->s;
}

int
rk_reset (void *vstate, size_t dim)
{
  rk_state_t *state = (rk_state_t *) vstate;

  cblas_dcopy (dim, 0, 0, state->y_onestep, 1);
  cblas_dcopy (dim, 0, 0, state->y0_orig, 1);
  cblas_dcopy (dim, 0, 0, state->y0, 1);
  cblas_dcopy (dim * state->s, 0, 0, state->fY, 1);
  cblas_dcopy (dim * state->s, 0, 0, state->fY_tmp, 1);

  return GSL_SUCCESS;
}
