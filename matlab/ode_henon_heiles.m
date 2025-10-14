function [f,H] = ode_henon_heiles(y)
p = y(:,1:2);
x = y(:,3:4);
dpdt(:,1) = -x(:,1) - 2 * x(:,1) .* x(:,2);
dpdt(:,2) = -x(:,2) - x(:,1) .^ 2 + x(:,2) .^ 2;
f = [dpdt,p];
if (nargout >= 2)
  H = sum(y(:,1:4) .^ 2,2) / 2 + x(:,1:2) .^ 2 * [1;-1/3] .* x(:,2);
end
