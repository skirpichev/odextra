#ifndef _ODEIV_UTILS_H_
#define _ODEIV_UTILS_H_

#include <gsl/gsl_odeiv.h>

typedef struct
{
  double *y0;
  double *y_onestep;
  double *y0_orig;
  size_t s; double *a, *b, *c;
  double *fY, *fY_tmp;
} rk_state_t;

rk_state_t*
rk_init (size_t dim, const double a[], const double b[],
         const double c[], const size_t s);

int
rk_apply (void *vstate, size_t dim, double t, double h, double y[],
          double yerr[], const double dydt_in[],
          double dydt_out[], const gsl_odeiv_system * sys);

unsigned int
rk_order (void *vstate);

void
rk_free (void *vstate);

int
rk_reset (void *vstate, size_t dim);

#endif /* _ODEIV_UTILS_H_ */
