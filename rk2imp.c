/* rk2imp.c */

/* Runge-Kutta 2, Gaussian implicit.  Also known as the implicit
   midpoint rule. */

#include <config.h>

#include "odeiv_util.h"

static void *
rk2imp_alloc (size_t dim)
{
  const double a[1] = { 0.5 };
  const double b[1] = { 1 };
  const double c[1] = { 0.5 };

  return rk_init (dim, a, b, c, 1);
}

static const gsl_odeiv_step_type rk2imp_type = {
  "rk2imp",      /* name */
  1,             /* can use dydt_in */
  1,             /* gives exact dydt_out */
  &rk2imp_alloc,
  &rk_apply,
  &rk_reset,
  &rk_order,
  &rk_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_midpoint = &rk2imp_type;
