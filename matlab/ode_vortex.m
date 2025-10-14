function [f,H,L,P,Q] = ode_vortex(y,gamma=[1,1,1,1])
z = p = y(:,1:4);
x = y(:,5:8);

dpdt(:,1:4) = zeros(rows(y),4);
dxdt(:,1:4) = zeros(rows(y),4);

for i=1:4
  z(:,i) = p(:,i) / gamma(i);
end

for i=1:4
  for j=1:4
    if i != j
      r = (x(:,i) - x(:,j)) .^2 + (z(:,i) - z(:,j)) .^2;
      dpdt(:,i) += gamma(i) * gamma(j) * (x(:,i) - x(:,j)) ./ r / 2;
      dxdt(:,i) -= gamma(j) * (z(:,i) - z(:,j)) ./ r / 2;
    end
  end
end

f = [dpdt,dxdt];
if nargout >= 2
  L = P = Q = H = zeros(rows(y),1);
  for i=1:4
    P += p(:,i);
    Q += gamma(i) * x(:,i);
    L += gamma(i) * (x(:,i) .^2 + z(:,i) .^2);
    for j=1:4
      if i != j
        r = sqrt((x(:,i) - x(:,j)) .^2 + (z(:,i) - z(:,j)) .^2);
        H -= gamma(i) * gamma(j) * log(r) / 4;
      end
    end
  end
end
