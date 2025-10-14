function [f,H] = ode_pendulum(y)
p = y(:,1);
x = y(:,2);
f = [-sin(x),p];
if (nargout == 2)
  H = p .^ 2 / 2 - cos(x);
end
