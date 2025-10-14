% From Hairer, p.242 (symplectic euler + orthogonal projection).

function demo_spherical_pendulum(dt=80)
y0 = [0.06,0,0,0,sin(0.1),-cos(0.1)];
te = linspace(0,dt,4e2);
tr = te;
ye = ceulers({@dH,@dG,@G},te,y0);
yr = rattle({@dH,@dG,@G},tr,y0);

subplot(2,1,1);
plot(te,G(ye), '-;Symplectic Euler with constraints;', ...
     tr,G(yr), '-;Rattle;');
title('Constraint');

subplot(2,1,2);
plot(te,H(ye) - H(y0), '-;Symplectic Euler with constraints;', ...
     tr,H(yr) - H(y0), '-;Rattle;');
title('Error in Hamiltonian');

figure;
plot3(ye(:,4), ye(:,5), ye(:,6), '-;Symplectic Euler with constraints;', ...
      yr(:,4), yr(:,5), yr(:,6), '-;Rattle;')

function dydt = dH(y)
dydt(:,4:6) = y(:,1:3);
dydt(:,1:2) = 0;
dydt(:,3) = -1;


function r = dG(y)
r = 2 * y(:,4:6);


function r = G(y)
r = sum(y(:,4:6).^2,2) - 1;


function r = H(y)
r = sum(y(:,1:3).^2,2) / 2 + y(:,6);
