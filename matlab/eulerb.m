% Implicit (backward) Euler integrator.

function y = eulerb(f,t,y0)
[h,y] = solver_init(t,y0);
abc = lobatto3b_abc(1);
y = solver_loop(h,y,@(h,y) rkstep(f,y,h,abc));
