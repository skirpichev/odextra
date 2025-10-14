function [h,y] = solver_init(tspan,init)
dim = length(init);
N = length(tspan);
h = diff(tspan);
y = [init;zeros(N-1,dim)];
