function [f,H] = ode_harmonic_oscillator(y,m = 1,k = 1)
p = y(:,1);
x = y(:,2);
f = [-k * x,p / m];
if (nargout >= 2)
  H = (p .^ 2 / m + k * x .^ 2) / 2;
end
