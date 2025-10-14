% Test for Henon-Heiles problem.

function demo_henon_heiles()
tspan = [0,15];
y0 = [1/6,0.1,0.1,0.1];
m = {'euler','eulerb','eulert','midpoint','eulers','verlet','gauss'};
printf('Demo: Henon-Heiles\n');
for i = 1:length(m)
  printf('solver: %s\n',m{i});
  switch m{i}
    case {'euler','eulerb'}
      N = 1e3;
    otherwise
      N = 1e2;
  end
  t{i} = linspace(tspan(1),tspan(2),N);
  y{i} = feval(m{i},@ode_henon_heiles,t{i},y0);
  [unused,H{i}] = ode_henon_heiles(y{i});
end
[unused,He] = ode_henon_heiles(y0);
figure;
semilogy(t{1},abs(H{1}-He)+eps, '-;Euler;', ...
         t{2},abs(H{2}-He)+eps, '-;Backward Euler;', ...
         t{3},abs(H{3}-He)+eps, '-;Trapesoidal rule;', ...
         t{4},abs(H{4}-He)+eps, '-;Implicit midpoint rule;', ...
         t{5},abs(H{5}-He)+eps, '-;Symplectic Euler;', ...
         t{6},abs(H{6}-He)+eps, '-;Verlet;', ...
         t{7},abs(H{7}-He)+eps, '-;Gauss;');
legend('location','southeast');
title('Henon-Heiles problem, hamiltonian conservation');
xlabel('t');
ylabel('|H-H(0)|');
