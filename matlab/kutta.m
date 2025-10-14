% Explicit classical Kutta method of order 4.

function y = kutta(f,t,y0)
[h,y] = solver_init(t,y0);
abc = rkinit_abc([0,0,0,0;1,0,0,0;0,1,0,0;0,0,2,0]/2,[1,2,2,1]/6);
y = solver_loop(h,y,@(h,y) explicit_rkstep(f,y,h,abc));
