% Verlet integrator.

function y = verlet(f,t,y0)
[h,y] = solver_init(t,y0);
abc = lobatto3a_abc(2);
ABC = lobatto3b_abc(2);
y = solver_loop(h,y,@(h,y) partitioned_rkstep(f,y,h,abc,ABC));

%!assert(verlet(@ode_pendulum,linspace(0,4 * pi),[0,2])(end,:),[0.0379597731890586,-1.9992078034945056],1e-2);
