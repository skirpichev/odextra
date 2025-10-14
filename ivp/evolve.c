/* evolve.c */

#include <config.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

gsl_odeiv_evolve *
gsl_odeiv_evolve_alloc (size_t dim)
{
  gsl_odeiv_evolve * e = (gsl_odeiv_evolve *) malloc (sizeof (gsl_odeiv_evolve));

  if (e == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for evolve struct", GSL_ENOMEM);
    }

  e->y0 = (double *) malloc (dim * sizeof (double));

  if (e->y0 == 0)
    {
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
    }

  e->yerr = (double *) malloc (dim * sizeof (double));

  if (e->yerr == 0)
    {
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for yerr", GSL_ENOMEM);
    }

  e->dydt_in = (double *) malloc (dim * sizeof (double));

  if (e->dydt_in == 0)
    {
      free (e->yerr);
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for dydt_in", GSL_ENOMEM);
    }

  e->dydt_out = (double *) malloc (dim * sizeof (double));

  if (e->dydt_out == 0)
    {
      free (e->dydt_in);
      free (e->yerr);
      free (e->y0);
      free (e);
      GSL_ERROR_NULL ("failed to allocate space for dydt_out", GSL_ENOMEM);
    }

  e->dimension = dim;
  e->count = 0;
  e->failed_steps = 0;
  e->last_step = 0.0;

  return e;
}

int
gsl_odeiv_evolve_reset (gsl_odeiv_evolve * e)
{
  e->count = 0;
  e->failed_steps = 0;
  e->last_step = 0.0;
  return GSL_SUCCESS;
}

void
gsl_odeiv_evolve_free (gsl_odeiv_evolve * e)
{
  RETURN_IF_NULL (e);
  free (e->dydt_out);
  free (e->dydt_in);
  free (e->yerr);
  free (e->y0);
  free (e);
}

/* Evolution framework method.
 *
 * Uses an adaptive step control object
 */
int
gsl_odeiv_evolve_apply (gsl_odeiv_evolve * e, gsl_odeiv_control * con, gsl_odeiv_step * step, const gsl_odeiv_system * dydt, double *t, double t1, double *h, double y[])
{
  const double t0 = *t;
  double h0 = *h;
  int step_status;
  int final_step = 0;
  double dt = t1 - t0; /* remaining time, possibly less than h */

  if (e->dimension != step->dimension)
    {
      GSL_ERROR ("step dimension must match evolution size", GSL_EINVAL);
    }

  if ((dt < 0.0 && h0 > 0.0) || (dt > 0.0 && h0 < 0.0))
    {
      GSL_ERROR ("step direction must match interval direction", GSL_EINVAL);
    }

  /* Save y in case of failure in a step.  No need to copy if we
     cannot control the step size. */

  if (con != NULL)
    {
      DBL_MEMCPY (e->y0, y, e->dimension);
    }

  /* Calculate initial dydt once or reuse previous value if the method
     can benefit. */

  if (step->type->can_use_dydt_in)
    {
      if (e->count == 0)
        {
          int status = GSL_ODEIV_FN_EVAL (dydt, t0, y, e->dydt_in);

          if (status)
            {
              return status;
            }
        }
      else
        {
          DBL_MEMCPY (e->dydt_in, e->dydt_out, e->dimension);
        }
    }

 try_step:

  if ((dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt))
    {
      h0 = dt;
      final_step = 1;
    }
  else
    {
      final_step = 0;
    }

  if (step->type->can_use_dydt_in)
    {
      step_status = gsl_odeiv_step_apply (step, t0, h0, y, e->yerr, e->dydt_in, e->dydt_out, dydt);
    }
  else
    {
      step_status = gsl_odeiv_step_apply (step, t0, h0, y, e->yerr, NULL, e->dydt_out, dydt);
    }

  /* Return if stepper indicates a pointer or user function failure */

  if (step_status == GSL_EFAULT || step_status == GSL_EBADFUNC)
    {
      return step_status;
    }

  /* Check for stepper internal failure */

  if (step_status != GSL_SUCCESS)
    {
      /* Stepper was unable to calculate step.  Try decreasing step size. */

      double h_old = h0;

      h0 *= 0.5;

      /* Check that an actual decrease in h0 occured and the
         suggested h0 will change the time by at least 1 ulp */

      {
        double t_curr = GSL_COERCE_DBL (*t);
        double t_next = GSL_COERCE_DBL ((*t) + h0);

        if (fabs (h0) < fabs (h_old) && t_next != t_curr)
          {
            /* Step was decreased.  Undo step, and try again with new h0. */

            DBL_MEMCPY (y, e->y0, dydt->dimension);
            e->failed_steps++;
            goto try_step;
          }
        else
          {
            *h = h0; /* notify user of step-size which caused the failure */
            *t = t0; /* restore original t value */
            return step_status;
          }
      }
    }

  e->count++;
  e->last_step = h0;

  if (final_step)
    {
      *t = t1;
    }
  else
    {
      *t = t0 + h0;
    }

  if (con != NULL)
    {
      /* Check error and attempt to adjust the step. */

      double h_old = h0;

      const int hadjust_status = gsl_odeiv_control_hadjust (con, step, y, e->yerr, e->dydt_out, &h0);

      if (hadjust_status == GSL_ODEIV_HADJ_DEC)
        {
          /* Check that the reported status is correct (i.e. an actual
             decrease in h0 occured) and the suggested h0 will change
             the time by at least 1 ulp */

          double t_curr = GSL_COERCE_DBL (*t);
          double t_next = GSL_COERCE_DBL ((*t) + h0);

          if (fabs (h0) < fabs (h_old) && t_next != t_curr)
            {
              /* Step was decreased.  Undo step, and try again with new h0. */

              DBL_MEMCPY (y, e->y0, dydt->dimension);
              e->failed_steps++;
              goto try_step;
            }
          else
            {
              /* Can not obtain required error tolerance, and can not
                 decrease step-size any further, so give up and return
                 GSL_FAILURE. */

              *h = h0; /* notify user of step-size which caused the failure */
              return GSL_FAILURE;
            }
        }
    }

  /* Suggest step size for next time-step.  Change of step size is not
     suggested in the final step, because that step can be very small
     compared to previous step, to reach t1. */

  if (final_step == 0)
    {
      *h = h0;
    }

  return step_status;
}
