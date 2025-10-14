function y = ceulers(f,t,y0)
[h,y] = solver_init(t,y0);

s = columns(y) / 2;
p = 1 : s;
x = p + s;

dH = f{1};
G = f{2};
g = f{3};

m = columns(g(y0));

for n = 2:rows(y)
  % Implicit RK step in p:
  y(n,p) = solver_solve(@(z) [z(p) - y(n-1,p) - ...
                              h(n-1) * dH([z(p),y(n-1,x)])(p) - ...
                              h(n-1) * G(y(n-1,:)) * z(s + 1:end), ...
                              g(y(n-1,:) + h(n-1) * dH([z(p),y(n-1,x)]))], ...
                        [y(n-1,p),zeros(1,m)])(p);

  % Explicit RK step in x:
  y(n,x) = y(n-1,x) + h(n-1) * dH([y(n,p),y(n-1,x)])(x);

  % Projection step:
  y(n,p) = solver_solve(@(z) [z(p) - y(n,p) - ...
                              h(n-1) * G(y(n,:)) * z(s + 1:end), ...
                              dH([z(p),y(n,x)])(x) * G(y(n,:))'], ...
                        [y(n,p),zeros(1,m)])(p);
end
