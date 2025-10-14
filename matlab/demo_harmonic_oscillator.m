% Test with harmonic oscillator.

function demo_harmonic_oscillator()
tspan = [0,20];
y0 = [1,0];
m = {'euler','eulerb','eulert','midpoint','eulers','verlet','gauss'};
printf('Demo: Harmonic oscillator\n');
for i = 1:length(m)
  printf('solver: %s\n',m{i});
  switch m{i}
    case {'euler','eulerb'}
      N = 1e3;
    otherwise
      N = 1e2;
  end
  t{i} = linspace(tspan(1),tspan(2),N);
  y{i} = feval(m{i},@ode_harmonic_oscillator,t{i},y0);
end
ye = [y0(1) * cos(tspan') - y0(2) * sin(tspan'), ...
      y0(1) * sin(tspan') + y0(2) * cos(tspan')];
figure;
plot(t{1},y{1}(:,1), '-;Euler;', ...
     t{2},y{2}(:,1), '-;Backward Euler;', ...
     t{3},y{3}(:,1), '-;Trapesoidal rule;', ...
     t{4},y{4}(:,1), '-;Implicit midpoint rule;', ...
     t{5},y{5}(:,1), '-;Symplectic Euler;', ...
     t{6},y{6}(:,1), '-;Verlet;', ...
     t{7},y{7}(:,1), '-;Gauss;', ...
     tspan',ye(:,1), 'x;Exact solution;');
title('Harmonic Oscillator');
xlabel('t');
ylabel('x');
