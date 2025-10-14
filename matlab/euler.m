% Explicit Euler integrator.

function y = euler(f,t,y0)
[h,y] = solver_init(t,y0);
abc = lobatto3a_abc(1);
y = solver_loop(h,y,@(h,y) explicit_rkstep(f,y,h,abc));
