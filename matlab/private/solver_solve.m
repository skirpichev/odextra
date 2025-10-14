function y = solver_solve(f,y0)
tol = eps(class(y0));
for n=1:400
  y = y0 - f(y0);
  if norm(y - y0,Inf) < tol
    break;
  end
  y0 = y;
end
