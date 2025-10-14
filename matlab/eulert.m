% Implicit (trapezoid rule) Euler integrator (Heun's method).

function y = eulert(f,t,y0)
[h,y] = solver_init(t,y0);
abc = lobatto3a_abc(2);
y = solver_loop(h,y,@(h,y) rkstep(f,y,h,abc));
