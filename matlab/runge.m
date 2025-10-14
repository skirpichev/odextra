% Explicit Runge method of order 2.

function y = runge(f,t,y0)
[h,y] = solver_init(t,y0);
abc = rkinit_abc([0,0;0.5,0],[0,1]);
y = solver_loop(h,y,@(h,y) explicit_rkstep(f,y,h,abc));
