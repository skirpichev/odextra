% Compute coefficients for Radau IIA method.

function abc = radau2a_abc(s)
abc = collocation_abc(legendre_coeff(s) - [0,legendre_coeff(s - 1)]);

%!assert(radau2a_abc(1),rkinit_abc(1,1),eps);
%!assert(radau2a_abc(2),rkinit_abc([5/3,-1/3;3,1]/4,[3,1]/4),eps);
%!assert(radau2a_abc(3),rkinit_abc([(88-7*sqrt(6))/360,(296-169*sqrt(6))/1800,(-2+3*sqrt(6))/225;(296+169*sqrt(6))/1800,(88+7*sqrt(6))/360,(-2-3*sqrt(6))/225;(16-sqrt(6))/36,(16+sqrt(6))/36,1/9],[(16-sqrt(6))/36,(16+sqrt(6))/36,1/9]),10*eps);
