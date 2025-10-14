#ifndef _ODE_TEST_H_
#define _ODE_TEST_H_

typedef int integral_func (double t, const double y[],
                           double integrals[], void *params);

typedef int exact_sol_func (const double t, double y[],
                            const double y0[], void *params);

GSL_VAR const gsl_odeiv_system *rhs_func_henon;
GSL_VAR integral_func int_henon;

GSL_VAR const gsl_odeiv_system *rhs_func_pendulum;
GSL_VAR integral_func int_pendulum;
GSL_VAR exact_sol_func exact_pendulum;

GSL_VAR const gsl_odeiv_system *rhs_func_kepler;
GSL_VAR integral_func int_kepler;
GSL_VAR exact_sol_func exact_kepler;

GSL_VAR const gsl_odeiv_system *rhs_func_linear;
GSL_VAR integral_func int_linear;
GSL_VAR exact_sol_func exact_linear;

GSL_VAR const gsl_odeiv_system *rhs_func_toda;
GSL_VAR integral_func int_toda;

GSL_VAR const gsl_odeiv_system *rhs_func_hill;
GSL_VAR integral_func int_hill;

GSL_VAR const gsl_odeiv_system *rhs_func_vortex;
GSL_VAR integral_func int_vortex;

GSL_VAR const gsl_odeiv_system *rhs_func_fermi;
GSL_VAR integral_func int_fermi;

#endif /* _ODE_TEST_H_ */
