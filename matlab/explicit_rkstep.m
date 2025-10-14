function dy = explicit_rkstep(f,y,h,abc)
k = ones(abc.s,1) * y;
ha = h * abc.a;
for i = 1:abc.s
  k(i,:) = f(y + ha(i,:) * k);
end
dy = h * abc.b * k;
