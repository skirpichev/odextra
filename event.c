/* event.c */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>

#include "odextra.h"
#include "tests/odextra_test.h"
#include "odeiv_util.h"

static gsl_rng *r;

static int
init (const double E, double y[])
{
  do
    {
      y[1] = 2.0 * gsl_rng_uniform (r) - 1;
      y[2] = 2.0 * gsl_rng_uniform (r) - 1;
      y[3] = 2.0 * gsl_rng_uniform (r) - 1;

      y[0] = 2.0 * E - y[1] * y[1] - y[2] * y[2] - y[3] * y[3] -
        2.0 * y[3] * (y[2] * y[2] - y[3] * y[3] / 3.0);
    }
  while (y[0] < 0);

  y[0] = GSL_SIGN (gsl_rng_uniform (r) - 0.5) * sqrt (y[0]);

  return GSL_SUCCESS;
}

static int
event_henon (const double t, const double y[], double* e, void *p)
{
  *e = y[2];

  return GSL_SUCCESS;
}

int
main (int argc, char *argv[])
{
  if (argc < 4)
    {
      return 1;
    }

  const gsl_odeiv_step_type *T = gsl_odeiv_step_midpoint;

  gsl_rng_env_setup();

  r = gsl_rng_alloc (gsl_rng_default);

  const double E = atof (argv[1]), t1 = atof (argv[2]);
  const size_t n = atoi (argv[3]);

  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_control *c = gsl_odeiv_control_fixed_new ();
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (4);

  for (size_t i = 0; i < n; i++)
    {
      double t = 0.0;
      double h = 1e-3;
      double y[4];

      init (E, y);

      gsl_odeiv_step_reset (s);
      gsl_odeiv_evolve_reset (e);

      while (t < t1)
        {
          double event[2];

          event_henon (t, y, &event[0], 0);

          int status;

          status = gsl_odeiv_evolve_apply (e, c, s,
                                           rhs_func_henon,
                                           &t, t1, &h, y);

          if (status != GSL_SUCCESS)
            {
              fprintf (stderr, "break");
              break;
            }

          /*
            dir determines whether events are detected only when the
            fold function is increasing (+1) or decreasing (-1).  In
            this case, we set dir = 0 because we want to stop in
            either case.
          */

          int dir = 1;

          event_henon (t, y, &event[1], 0);

          if (event[0] * event[1] <= 0 && dir * event[1] <= 0)
            {
              fprintf (stdout, "%+.7e\t%+.7e\t%+.7e\t%+.7e\t%+.7e\n",
                       t, y[0], y[1], y[2], y[3]);
            }
        }
    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  gsl_rng_free (r);

  return 0;
}

/*
  Local Variables:
  compile-command: "c99 -I. -L ./.libs/ -lgsl -lgslcblas \
  event.c rk2imp.c henon.c cfxd.c rkstep.c -o event"
  End:
*/
