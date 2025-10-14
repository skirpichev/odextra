% Test matrix a, if it's an explicit RK method.

function r = isexplrk(a)
r = issquare(a) && all(a == tril(a,-1));

%!assert(isexplrk(0),true);
%!assert(isexplrk(1),false);
%!assert(isexplrk(tril(rand(4),-1)),true);
