/* control.c */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

gsl_odeiv_control *
gsl_odeiv_control_alloc (const gsl_odeiv_control_type * T)
{
  gsl_odeiv_control * c = (gsl_odeiv_control *) malloc (sizeof (gsl_odeiv_control));

  if(c == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for control struct", GSL_ENOMEM);
    }

  c->type = T;
  c->state = c->type->alloc ();

  if (c->state == 0)
    {
      free (c); /* exception in constructor, avoid memory leak */

      GSL_ERROR_NULL ("failed to allocate space for control state", GSL_ENOMEM);
    }

  return c;
}

int
gsl_odeiv_control_init (gsl_odeiv_control * c, double eps_abs, double eps_rel, double a_y, double a_dydt)
{
  return c->type->init (c->state, eps_abs, eps_rel, a_y, a_dydt);
}

void
gsl_odeiv_control_free (gsl_odeiv_control * c)
{
  RETURN_IF_NULL (c);
  c->type->free (c->state);
  free (c);
}

const char *
gsl_odeiv_control_name (const gsl_odeiv_control * c)
{
  return c->type->name;
}

int
gsl_odeiv_control_hadjust (gsl_odeiv_control * c, gsl_odeiv_step * s, const double y[], const double yerr[], const double dydt[], double * h)
{
  return c->type->hadjust (c->state, s->dimension, s->type->order (s->state), y, yerr, dydt, h);
}
