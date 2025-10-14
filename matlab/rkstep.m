function dy = rkstep(f,y,h,abc)
k = ones(abc.s,1) * y;
ha = h * abc.a;
dy = h * abc.b * solver_solve(@(_) _ - f(k + ha * _),k);
