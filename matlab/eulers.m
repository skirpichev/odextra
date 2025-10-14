% Symplectic Euler integrator.

function y = eulers(f,t,y0)
[h,y] = solver_init(t,y0);
abc = lobatto3a_abc(1);
ABC = lobatto3b_abc(1);
y = solver_loop(h,y,@(h,y) partitioned_rkstep(f,y,h,abc,ABC));
