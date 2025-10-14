% Compute coefficients for Lobatto IIIA method.

function abc = lobatto3a_abc(s)
abc = collocation_abc([legendre_coeff(s - 1) ./ (s:-1:1),0]);

%!assert(lobatto3a_abc(1),rkinit_abc(0,1),eps);
%!assert(lobatto3a_abc(2),rkinit_abc([0,0;0.5,0.5],[0.5,0.5]),eps);
%!assert(lobatto3a_abc(3),rkinit_abc([0,0,0;5/24,1/3,-1/24;1/6,2/3,1/6],[1/6,2/3,1/6]),eps);
%!assert(lobatto3b_abc(4),rkinit_abc([1/12,(-1-sqrt(5))/24,(-1+sqrt(5))/24,0;1/12,(25+sqrt(5))/120,(25-13*sqrt(5))/120,0;1/12,(25+13*sqrt(5))/120,(25-sqrt(5))/120,0;1/12,(11-sqrt(5))/24,(11+sqrt(5))/24,0],[1/12,5/12,5/12,1/12]),100*eps);
