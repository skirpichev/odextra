function dy = partitioned_rkstep(f,y,h,abc,ABC,p = columns(y) / 2,x = p)
k = ones(abc.s,1) * y;
ha = h * abc.a;
hA = h * ABC.a;
x = p + (1:x);
p = 1:p;
k = solver_solve(@(_) _ - f(k + [ha * _(:,p),hA * _(:,x)]),k);
dy = h * [abc.b * k(:,p),ABC.b * k(:,x)];
