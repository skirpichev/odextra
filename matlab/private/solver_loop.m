function y = solver_loop(h,y,dy)
c = zeros(size(y(1,:)));
for n = 2:rows(y)
  % Use compensated summation:
  c += dy(h(n - 1),y(n - 1,:));
  y(n,:) = y(n - 1,:) + c;
  c += y(n - 1,:) - y(n,:);
end
