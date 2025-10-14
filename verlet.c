/* verlet.c */

#include <config.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include "odextra.h"
#include "odeiv_util.h"

/* Verlet */

typedef struct
{
  double *f;
  double *f_half;
  double *f2;
  double *y0;
  double *y_tmp;
  double *y_onestep;
}
  verlet_state_t;

static void *
verlet_alloc (size_t dim)
{
  verlet_state_t *state = (verlet_state_t *) malloc (sizeof (verlet_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for state", GSL_ENOMEM);
    }

  state->f = (double *) malloc (dim * sizeof (double));

  if (state->f == 0)
    {
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for f", GSL_ENOMEM);
    }

  state->f_half = (double *) malloc (dim * sizeof (double));

  if (state->f_half == 0)
    {
      free (state->f);
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for f_half", GSL_ENOMEM);
    }

  state->f2 = (double *) malloc (dim * sizeof (double));

  if (state->f2 == 0)
    {
      free (state->f);
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for f_half", GSL_ENOMEM);
    }

  state->y0 = (double *) malloc (dim * sizeof (double));

  if (state->y0 == 0)
    {
      free (state->f2);
      free (state->f_half);
      free (state->f);
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
    }

  state->y_tmp = (double *) malloc (dim * sizeof (double));

  if (state->y_tmp == 0)
    {
      free (state->y0);
      free (state->f2);
      free (state->f_half);
      free (state->f);
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for y_tmp", GSL_ENOMEM);
    }

  state->y_onestep = (double *) malloc (dim * sizeof (double));

  if (state->y_onestep == 0)
    {
      free (state->y_tmp);
      free (state->y0);
      free (state->f2);
      free (state->f_half);
      free (state->f);
      free (state);

      GSL_ERROR_NULL ("failed to allocate space for y_onestep", GSL_ENOMEM);
    }

  return state;
}

static int
verlet_step (double *y, const verlet_state_t * state,
             const double h, const double t, const size_t dim,
             const gsl_odeiv_system * sys)
{
  /* Make a symplectic Verlet advance with step h */

  const size_t n = dim / 2;

  double *y0 = state->y_tmp;

  cblas_dcopy (dim, y, 1, y0, 1); /* Save original y */

  double *f = state->f;

  /* Solve p_{n+1/2} = p_{n} - h / 2 * dH/dx(p_{n+1/2},x_{n}) */

  const int maxiter = 1e3;
  int iter = 1;

  const double normtol = 10 * GSL_DBL_EPSILON;
  double norm;

  do
    {
      /* Method is explicit when H(p,x) = K(p) + U(x) */

      norm = 0.0; /* reset norm */

      for (size_t i = 0; i < n; i++) /* p[0] = y[0] ... */
        {
          double p0 = y[i];

          y[i] = y0[i] + (h / 2) * f[i];

          double dp = fabs (y[i] - p0);

          norm = (dp > norm) ? dp : norm;
        }

      int s = GSL_ODEIV_FN_EVAL (sys, t, y, f);

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      iter++;

      if (iter > maxiter)
        {
          GSL_ERROR ("reached iteration limit for "
                     "fixed point equation solver", GSL_EMAXITER);
        }
    }
  while (norm > normtol);

  double *f2 = state->f2;

  cblas_dcopy (dim, f, 1, f2, 1);

  /* Solve x_{n+1} = x_{n} + h / 2 * (dH/dp(p_{n+1/2},x_{n}) +
     dH/dp(p_{n+1/2},x_{n+1})) */

  iter = 1;

  do
    {
      /* Method is explicit when H(p,x) = K(p) + U(x) */

      norm = 0.0; /* reset norm */

      for (size_t i = n; i < dim; i++) /* x[0] = y[n] ... */
        {

          double x0 = y[i];

          y[i] = y0[i] + (h / 2) * (f[i] + f2[i]);

          double dx = fabs (y[i] - x0);

          norm = (dx > norm) ? dx : norm;
        }

      int s = GSL_ODEIV_FN_EVAL (sys, t + h, y, f2);

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      iter++;

      if (iter > maxiter)
        {
          GSL_ERROR ("reached iteration limit for "
                     "fixed point equation solver", GSL_EMAXITER);
        }
    }
  while (norm > normtol);

  /* Explicit step: Advance in p with step size h */

  for (size_t i = 0; i < n; i++) /* p[0] = y[0] ... */
    {
      y[i] = y0[i] + (h / 2) * (f[i] + f2[i]);
    }

  return GSL_SUCCESS;
}

static int
verlet_apply (void *vstate, size_t dim, double t, double h,
              double y[], double yerr[], const double dydt_in[],
              double dydt_out[], const gsl_odeiv_system * sys)
{
  verlet_state_t *state = (verlet_state_t *) vstate;

  double *y0 = state->y0;

  cblas_dcopy (dim, y, 1, y0, 1); /* save original y for possible failures */

  /* Calculate f = f(t, y) */

  double *f = state->f;

  if (dydt_in != NULL)
    {
      cblas_dcopy (dim, dydt_in, 1, f, 1);
    }
  else
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t, y, f);

      if (s != GSL_SUCCESS)
        {
          return s;
        }
    }

  double *f_half = state->f_half;

  cblas_dcopy (dim, f, 1, f_half, 1); /* Save f values for half-step stage */

  /* First traverse h with one step (save to y_onestep) */

  double *y_onestep = state->y_onestep;

  cblas_dcopy (dim, y, 1, y_onestep, 1);

  {
    int s = verlet_step (y_onestep, state, h, t, dim, sys);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* Advance half step for error estimation */

  cblas_dcopy (dim, f_half, 1, f, 1); /* Restore previous f values */

  /* Half-step one */

  {
    int s = verlet_step (y, state, h / 2.0, t, dim, sys);

    if (s != GSL_SUCCESS)
      {
        cblas_dcopy (dim, y0, 1, y, 1); /* restore original y */

        return s;
      }
  }

  /* Half-step two */

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h / 2.0, y, f);

    if (s != GSL_SUCCESS)
      {
        return s;
      }

    s = verlet_step (y, state, h / 2.0, t + h / 2.0, dim, sys);

    if (s != GSL_SUCCESS)
      {
        cblas_dcopy (dim, y0, 1, y, 1); /* restore original y */

        return s;
      }
  }

  /* Derivatives at output */

  if (dydt_out != NULL)
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t + h, y, dydt_out);

      if (s != GSL_SUCCESS)
        {
          cblas_dcopy (dim, y0, 1, y, 1); /* restore original y */

          return s;
        }
    }

  /* Error estimation

     yerr = C * 0.5 * | y(onestep) - y(twosteps) | / (2^order - 1)

     Constant C is approximately 8.0 to ensure 90% of samples lie
     within the error (assuming a gaussian distribution with prior
     p(sigma)=1/sigma.)
  */

  for (size_t i = 0; i < dim; i++)
    {
      yerr[i] = 4.0 * (y[i] - y_onestep[i]) / 3.0;
    }

  return GSL_SUCCESS;
}

static int
verlet_reset (void *vstate, size_t dim)
{
  verlet_state_t *state = (verlet_state_t *) vstate;

  cblas_dcopy (dim, 0, 0, state->f, 1);
  cblas_dcopy (dim, 0, 0, state->f_half, 1);
  cblas_dcopy (dim, 0, 0, state->f2, 1);
  cblas_dcopy (dim, 0, 0, state->y0, 1);
  cblas_dcopy (dim, 0, 0, state->y_tmp, 1);
  cblas_dcopy (dim, 0, 0, state->y_onestep, 1);

  return GSL_SUCCESS;
}

static unsigned int
verlet_order (void *vstate)
{
  (void) vstate;

  return 2;
}

static void
verlet_free (void *vstate)
{
  verlet_state_t *state = (verlet_state_t *) vstate;

  free (state->f);
  free (state->f_half);
  free (state->f2);
  free (state->y0);
  free (state->y_tmp);
  free (state->y_onestep);
  free (state);
}

static const gsl_odeiv_step_type verlet_type = {
  "verlet",      /* name */
  1,             /* can use dydt_in */
  1,             /* gives exact dydt_out */
  &verlet_alloc, /* state allocation */
  &verlet_apply, /* advance one step */
  &verlet_reset, /* reset to well-known state */
  &verlet_order, /* return the order of the stepping function */
  &verlet_free   /* free state */
};

const gsl_odeiv_step_type *gsl_odeiv_step_verlet = &verlet_type;
