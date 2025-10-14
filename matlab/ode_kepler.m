function [f,H,L] = ode_kepler(y)
p = y(:,1:2);
x = y(:,3:4);
r = sqrt(x(:,1) .^ 2 + x(:,2) .^ 2);
dpdt(:,1) = -x(:,1) ./ r .^ 3;
dpdt(:,2) = -x(:,2) ./ r .^ 3;
f = [dpdt,p];
if (nargout >= 2)
  H = (p(:,1) .^ 2 + p(:,2) .^ 2) / 2 - 1 ./ r;
end
if (nargout >= 3)
  L = x(:,1) .* p(:,2) - x(:,2) .* p(:,1);
end
