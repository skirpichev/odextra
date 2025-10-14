function y = solver_fsolve(f,y)
options = optimset('TolX',eps,'TolFun',eps);
[y,unused,status,info] = fsolve(f,y,options);
if (status != 1 && status != -3) || (status == -3 && info.iterations > 1)
  printf('status %d in %d iterations\n',status,info.iterations);
end
