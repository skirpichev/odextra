% Implicit (midpoint rule) Euler integrator.

function y = midpoint(f,t,y0)
[h,y] = solver_init(t,y0);
abc = gauss_abc(1);
y = solver_loop(h,y,@(h,y) rkstep(f,y,h,abc));
