/* step.c */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

gsl_odeiv_step *
gsl_odeiv_step_alloc (const gsl_odeiv_step_type * T, size_t dim)
{
  gsl_odeiv_step *s = (gsl_odeiv_step *) malloc (sizeof (gsl_odeiv_step));

  if (s == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for ode struct", GSL_ENOMEM);
    }

  s->type = T;
  s->dimension = dim;

  s->state = s->type->alloc (dim);

  if (s->state == 0)
    {
      free (s); /* exception in constructor, avoid memory leak */

      GSL_ERROR_NULL ("failed to allocate space for ode state", GSL_ENOMEM);
    }

  return s;
}

const char *
gsl_odeiv_step_name (const gsl_odeiv_step * s)
{
  return s->type->name;
}

unsigned int
gsl_odeiv_step_order (const gsl_odeiv_step * s)
{
  return s->type->order (s->state);
}

int
gsl_odeiv_step_apply (gsl_odeiv_step * s, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)
{
  return s->type->apply (s->state, s->dimension, t, h, y, yerr, dydt_in, dydt_out, dydt);
}

int
gsl_odeiv_step_reset (gsl_odeiv_step * s)
{
  return s->type->reset (s->state, s->dimension);
}

void
gsl_odeiv_step_free (gsl_odeiv_step * s)
{
  RETURN_IF_NULL (s);
  s->type->free (s->state);
  free (s);
}
