function [f,H,P,L] = ode_solar(y,m,G = 0.01720209895^2)
n = length(m);
p = y(:,1:(3 * n));
v = p;
x = y(:,(3 * n + 1):end);

dpdt = zeros(rows(y),3 * n);

for i=1:n
  v(:,(3 * i - 2):(3 * i)) = p(:,(3 * i - 2):(3 * i)) / m(i);
  for j=1:n
    if i != j
      r = x(:,(3 * i - 2):(3 * i)) - x(:,(3 * j - 2):(3 * j));
      dpdt(:,(3 * i - 2):(3 * i)) -= G * m(i) * m(j) * r ...
          ./ ((sum(r.^2,2).^(3/2)) * ones(1,3));
    end
  end
end

f = [dpdt,v];

if nargout >= 2
  H = zeros(rows(y),1);
  P = zeros(rows(y),3);
  L = zeros(rows(y),3);
  for i=1:n
    P += p(:,(3 * i - 2):(3 * i));
    L += [p(:,3 * i - 1) .* x(:,3 * i) - ...
          p(:,3 * i) .* x(:,3 * i - 1), ...
          p(:,3 * i) .* x(:,3 * i - 2) - ...
          p(:,3 * i - 2) .* x(:,3 * i), ...
          p(:,3 * i - 2) .* x(:,3 * i - 1) - ...
          p(:,3 * i - 1) .* x(:,3 * i - 2)];
    H += sum(p(:,(3 * i - 2):(3 * i)).^2,2) / m(i) / 2;
    for j=1:n
      if i > j
        r = x(:,(3 * i - 2):(3 * i)) - x(:,(3 * j - 2):(3 * j));
        H -= G * m(i) * m(j) ./ sqrt(sum(r.^2,2));
      end
    end
  end
end
