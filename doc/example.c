#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>

#include <odextra.h>

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

static const gsl_odeiv_system sys_henon = { rhs_henon, 0, 4, 0 };


int
main (int argc, char *argv[])
{
  (void) argc;
  (void) argv;

  const gsl_odeiv_step_type *T = gsl_odeiv_step_gauss;
  const double t1 = 10;

  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (4);

  double y[4] = {1/6,0.1,0.1,0.1}, t = 0.0, h = 1e-3;

  while (t < t1)
    {
      printf ("%.5e %+.5e %+.5e %+.5e %+.5e\n", t, y[0], y[1], y[2], y[3]);

      int status = gsl_odeiv_evolve_apply (e, NULL, s,
                                           &sys_henon,
                                           &t, t1, &h, y);

      if (status != GSL_SUCCESS)
        {
          fprintf (stderr,
                   "error, return value=%d\n",
                   status);
          break;
        }
    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_step_free (s);

  return 0;
}
