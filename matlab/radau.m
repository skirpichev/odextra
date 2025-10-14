% Radau IIA of order 5.

function y = radau(f,t,y0)
[h,y] = solver_init(t,y0);
abc = radau2a_abc(3);
y = solver_loop(h,y,@(h,y) rkstep(f,y,h,abc));

%!assert(radau(@ode_pendulum,linspace(0,4 * pi),[0,2])(end,:),[0.0379597731890586,-1.9992078034945056],1e-5);
