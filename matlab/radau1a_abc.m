% Compute coefficients for Radau IA method.

function abc = radau1a_abc(s)
abc = collocation_abc(legendre_coeff(s) + [0,legendre_coeff(s - 1)]);
V = vander(abc.c);
V = V(:,abc.s:-1:1);
B = abc.b' * ones(1,abc.s);
J = B' .* V';
F = diag(1 ./ (1:abc.s));
D = ones(abc.s) - (diag(abc.c) * V)';
abc.a = J^(-1) * (F * (B' .* D));

%!assert(radau1a_abc(1),rkinit_abc(1,1,0));
%!assert(radau1a_abc(2),rkinit_abc([1,-1;1,5/3]/4,[1,3]/4),eps);
%!assert(radau1a_abc(3),rkinit_abc([1/9,(-1-sqrt(6))/18,(-1+sqrt(6))/18;1/9,(88+7*sqrt(6))/360,(88-43*sqrt(6))/360;1/9,(88+43*sqrt(6))/360,(88-7*sqrt(6))/360],[1/9,(16+sqrt(6))/36,(16-sqrt(6))/36]),10*eps);
