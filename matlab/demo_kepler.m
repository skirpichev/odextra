% Test for Kepler problem.

function demo_kepler()
tspan = [0,12];
y0 = [0,1,1,0];
m = {'euler','eulerb','eulert','midpoint','eulers','verlet','gauss'};
printf('Demo: Kepler\n');
for i = 1:length(m)
  printf('solver: %s\n',m{i});
  switch m{i}
    case {'euler','eulerb'}
      N = 1e3;
    otherwise
      N = 1e2;
  end
  t{i} = linspace(tspan(1),tspan(2),N);
  y{i} = feval(m{i},@ode_kepler,t{i},y0);
  [f,H{i},M{i}] = ode_kepler(y{i});
end
[unused,He,Me] = ode_kepler(y0);
figure;
semilogy(t{1},abs(M{1}-Me)+eps, '-;Euler;', ...
         t{2},abs(M{2}-Me)+eps, '-;Backward Euler;', ...
         t{3},abs(M{3}-Me)+eps, '-;Trapesoidal rule;', ...
         t{4},abs(M{4}-Me)+eps, '-;Implicit midpoint rule;', ...
         t{5},abs(M{5}-Me)+eps, '-;Symplectic Euler;', ...
         t{6},abs(M{6}-Me)+eps, '-;Verlet;',...
         t{7},abs(M{7}-Me)+eps, '-;Gauss;');
legend('location','east');
title('Kepler problem, angular momentum conservation');
xlabel('t');
ylabel('|M-M(0)|');
figure;
semilogy(t{1},abs(H{1}-He)+eps, '-;Euler;', ...
         t{2},abs(H{2}-He)+eps, '-;Backward Euler;', ...
         t{3},abs(H{3}-He)+eps, '-;Trapesoidal rule;', ...
         t{4},abs(H{4}-He)+eps, '-;Implicit midpoint rule;', ...
         t{5},abs(H{5}-He)+eps, '-;Symplectic Euler;', ...
         t{6},abs(H{6}-He)+eps, '-;Verlet;', ...
         t{7},abs(H{7}-He)+eps, '-;Gauss;');
legend('location','southeast');
title('Kepler problem, hamiltonian conservation');
xlabel('t');
ylabel('|H-H(0)|');
