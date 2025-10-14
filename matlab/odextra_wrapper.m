function [t,y] = odextra_wrapper(f,t,y0,solver=@gauss)
fv = @(y) vecodefunc (@(y_) f (0, y_')', y);
y0=y0(:)';
if columns(t) == 2
  N = 1000;
  t = linspace(t(1),t(2),N);
end
y=solver(fv,t,y0);
t=t';
