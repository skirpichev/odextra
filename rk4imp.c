/* rk4imp.c */

/* Runge-Kutta 4, Gaussian implicit.  Also known as the
   Hammer-Hollingsworth method. */

#include <config.h>

#include <gsl/gsl_math.h>

#include "odeiv_util.h"

static void *
rk4imp_alloc (size_t dim)
{
  const double a[4] = { 0.25,                 0.25 - 0.5 / M_SQRT3,
                        0.25 + 0.5 / M_SQRT3, 0.25 };
  const double b[2] = { 0.5, 0.5 };
  const double c[2] = { 0.5 - 0.5 / M_SQRT3, 0.5 + 0.5 / M_SQRT3 };

  return rk_init (dim, a, b, c, 2);
}

static const gsl_odeiv_step_type rk4imp_type = {
  "rk4imp",      /* name */
  1,             /* can use dydt_in? */
  1,             /* gives exact dydt_out? */
  &rk4imp_alloc,
  &rk_apply,
  &rk_reset,
  &rk_order,
  &rk_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_gauss = &rk4imp_type;
