#ifndef _ODEXTRA_H_
#define _ODEXTRA_H_

#include <gsl/gsl_types.h>
#include <gsl/gsl_odeiv.h>

GSL_VAR const gsl_odeiv_step_type *gsl_odeiv_step_euler;
GSL_VAR const gsl_odeiv_step_type *gsl_odeiv_step_verlet;
GSL_VAR const gsl_odeiv_step_type *gsl_odeiv_step_midpoint;
GSL_VAR const gsl_odeiv_step_type *gsl_odeiv_step_gauss;

#endif /* _ODEXTRA_H_ */
