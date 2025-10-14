% Test with pendulum.

function demo_pendulum()
tspan = [0,15];
y0 = [1,0];
m = {'euler','eulerb','eulert','midpoint','eulers','verlet','gauss'};
printf('Demo: Pendulum\n');
for i = 1:length(m)
  printf('solver: %s\n',m{i});
  switch m{i}
    case {'euler','eulerb'}
      N = 1e3;
    otherwise
      N = 1e2;
  end
  t{i} = linspace(tspan(1),tspan(2),N);
  y{i} = feval(m{i},@ode_pendulum,t{i},y0);
  [unused,H{i}] = ode_pendulum(y{i});
end
[unused,He] = ode_pendulum(y0);
m = sqrt((1 + He) / 2);
[sn,cn,dn] = ellipj(tspan',m^2);
ye(:,1) = 2 * m * real(cn);
ye(:,2) = 2 * asin(m * real(sn));
figure;
plot(t{1},y{1}(:,1), '-;Euler;', ...
     t{2},y{2}(:,1), '-;Backward Euler;', ...
     t{3},y{3}(:,1), '-;Trapesoidal rule;', ...
     t{4},y{4}(:,1), '-;Implicit midpoint rule;', ...
     t{5},y{5}(:,1), '-;Symplectic Euler;', ...
     t{6},y{6}(:,1), '-;Verlet;', ...
     t{7},y{7}(:,1), '-;Gauss;',...
     tspan',ye(:,1), 'x;Exact solution;');
title('Pendulum');
xlabel('t');
ylabel('x');
figure;
semilogy(t{1},abs(H{1}-He)+eps, '-;Euler;', ...
         t{2},abs(H{2}-He)+eps, '-;Backward Euler;', ...
         t{3},abs(H{3}-He)+eps, '-;Trapesoidal rule;', ...
         t{4},abs(H{4}-He)+eps, '-;Implicit midpoint rule;', ...
         t{5},abs(H{5}-He)+eps, '-;Symplectic Euler;', ...
         t{6},abs(H{6}-He)+eps, '-;Verlet;', ...
         t{7},abs(H{7}-He)+eps, '-;Gauss;');
legend('location','southeast');
title('Pendulum, hamiltonian conservation');
xlabel('t');
ylabel('|H-H(0)|');
