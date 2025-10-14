/* test.c */

#include <config.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_test.h>

#include "odextra.h"
#include "odextra_test.h"
#include "odeiv_util.h"

static void
test_evolve_system (const gsl_odeiv_step_type * T,
                    const gsl_odeiv_system * sys, double t0, double t1,
                    double hstart, double y[], const char *desc)
{
  /* Tests system sys with stepper T.  Step length is controlled by
     error estimation from the stepper. */

  const size_t maxiter = 1e+5;

  int steps = 0;

  double t = t0;
  double h = hstart;

  gsl_odeiv_step *step = gsl_odeiv_step_alloc (T, sys->dimension);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (sys->dimension);

  double *y_orig = malloc (sys->dimension * sizeof (double));

  while (t < t1)
    {
      double t_orig = t;
      int s;

      cblas_dcopy (sys->dimension, y, 1, y_orig, 1);
      s = gsl_odeiv_evolve_apply (e, NULL, step, sys, &t, t1, &h, y);

      if (s != GSL_SUCCESS)
        {
          /* Check that t and y are unchanged */

          gsl_test_abs (t, t_orig, 0, "%s, t must be restored on failure",
                        gsl_odeiv_step_name (step));

          for (size_t i = 0; i < sys->dimension; i++)
            {
              gsl_test_abs (y[i], y_orig[i], 0,
                            "%s, y must be restored on failure",
                            gsl_odeiv_step_name (step), desc, i);
            }

          gsl_test (s, "%s evolve_apply returned %d",
                    gsl_odeiv_step_name (step), s);
          break;
        }

      if (steps > maxiter)
        {
          gsl_test (GSL_EFAILED,
                    "%s evolve_apply reached maxiter = %d",
                    gsl_odeiv_step_name (step), maxiter);
          break;
        }

      steps++;
    }

  free (y_orig);
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_step_free (step);
}

static void
test_evolve_pendulum (const gsl_odeiv_step_type * T, double h, double err)
{
  double y0[2] = { 1, 0 }, yf[2];
  const double t = 4 * M_PI;
  const size_t factor = 10;
  const char desc[] = "pendulum";

  exact_pendulum (t, yf, y0, NULL);

  test_evolve_system (T, rhs_func_pendulum, 0, t, h, y0, desc);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < rhs_func_pendulum->dimension; i++)
    {
      gsl_test_abs (y0[i], yf[i], factor * err,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}

static void
test_evolve_linear (const gsl_odeiv_step_type * T, double h, double err)
{
  double y0[4] = { 1, 0.5, 0, 0 }, yf[4];
  const double t = 4 * M_PI;
  const size_t factor = 10;
  const char desc[] = "linear";

  exact_linear (t, yf, y0, NULL);

  test_evolve_system (T, rhs_func_linear, 0, t, h, y0, desc);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < rhs_func_linear->dimension; i++)
    {
      gsl_test_abs (y0[i], yf[i], factor * err,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}

static void
test_evolve_kepler (const gsl_odeiv_step_type * T, double h, double err)
{
  double y0[4] = { 0, 1, 1, 0 }, yf[4];

  const double K = (y0[0] * y0[0] + y0[1] * y0[1]) / 2.0;
  const double U = -1.0 / sqrt (y0[2] * y0[2] + y0[3] * y0[3]);
  const double E = K + U;
  const double a = -1 / 2.0 / E;
  const double t = sqrt (gsl_pow_3 (a)) * 4 * M_PI;

  const size_t factor = 10;
  const char desc[] = "kepler";

  exact_kepler (t, yf, y0, NULL);

  test_evolve_system (T, rhs_func_kepler, 0, t, h, y0, desc);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < rhs_func_kepler->dimension; i++)
    {
      gsl_test_abs (y0[i], yf[i], factor * err,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}


static void
test_int_toda (const gsl_odeiv_step_type * T, double h)
{
  double y0[20];

  for (size_t i = 0; i < 10; i++)
    {
      y0[i] = 0;
      y0[i + 10] = sin (2 * M_PI * i / 8.0);
    }

  double I0[10], I[10];
  const double t = 4 * M_PI;
  const size_t factor = 10;
  const char desc[] = "toda";

  int_toda (0, y0, I0, 0);

  test_evolve_system (T, rhs_func_toda, 0, t, h, y0, desc);

  int_toda (t, y0, I, 0);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < 10; i++)
    {
      gsl_test_abs (I[i], I0[i], factor * h,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}

static void
test_int_vortex (const gsl_odeiv_step_type * T, double h)
{
  double y0[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  double I0[4], I[4];
  const double t = 4 * M_PI;
  const size_t factor = 10;
  const char desc[] = "vortex";
  double gamma[4] = {1, 2, 3, 4};

  int_vortex (0, y0, I0, gamma);

  test_evolve_system (T, rhs_func_vortex, 0, t, h, y0, desc);

  int_vortex (t, y0, I, gamma);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < 4; i++)
    {
      gsl_test_abs (I[i], I0[i], factor * h,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}

static void
test_int_fermi (const gsl_odeiv_step_type * T, double h)
{
  double y0[128];

  for (size_t i = 0; i < 64; i++)
    {
      y0[i] = 0;
      y0[i + 64] = sin (2 * M_PI * i / 63.0);
    }

  double I0[1], I[1];
  const double t = 4 * M_PI;
  const size_t factor = 10;
  const char desc[] = "fermi";
  double alpha = 0.25;

  int_fermi (0, y0, I0, &alpha);

  test_evolve_system (T, rhs_func_fermi, 0, t, h, y0, desc);

  int_fermi (t, y0, I, &alpha);

  /* err_target is target error of one step. Test if stepper has made
     larger error than (tolerance factor times) the number of steps
     times the err_target */

  for (size_t i = 0; i < 1; i++)
    {
      gsl_test_abs (I[i], I0[i], factor * h,
                    "%s %s evolve(%d)", T->name, desc, i);
    }
}

int
main (int argc, char *argv[])
{
  (void) argc, argv;

  struct stype
  {
    const gsl_odeiv_step_type *type;
    double h;
    double err;
  } s[5] = {
    { .type = gsl_odeiv_step_euler, .h = 4 * M_PI / 1e+4, .err = 4 * M_PI / 1e+4 },
    { .type = gsl_odeiv_step_verlet, .h = 4 * M_PI / 1e+3, .err = gsl_pow_2 (4 * M_PI / 1e+3) },
    { .type = gsl_odeiv_step_midpoint, .h = 4 * M_PI / 1e+3, .err = gsl_pow_2 (4 * M_PI / 1e+3) },
    { .type = gsl_odeiv_step_gauss, .h = 4 * M_PI / 1e+2, .err = gsl_pow_4 (4 * M_PI / 1e+2) },
    { .type = 0 }
  };

  gsl_ieee_env_setup ();

  for (size_t i = 1; s[i].type; i++)
    {
      test_evolve_pendulum (s[i].type, s[i].h, s[i].err);
      test_evolve_linear (s[i].type, s[i].h, s[i].err);
      test_evolve_kepler (s[i].type, s[i].h, s[i].err);

      test_int_toda (s[i].type, s[i].h);
      test_int_vortex (s[i].type, s[i].h / 11);
      test_int_fermi (s[i].type, s[i].h);
    }

  return gsl_test_summary ();
}
