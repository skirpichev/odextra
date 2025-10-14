% Compute coefficients for Lobatto IIIB method.

function abc = lobatto3b_abc(s)
abc = lobatto3a_abc(s);

abc.a = ones(s,1) * abc.b - ((1 ./ abc.b)' * abc.b) .* abc.a';
abc.c = sum(abc.a,2); % for s=1,2

%!assert(lobatto3b_abc(1),rkinit_abc(1,1),eps);
%!assert(lobatto3b_abc(2),rkinit_abc([0.5,0;0.5,0],[0.5,0.5]),eps);
%!assert(lobatto3b_abc(3),rkinit_abc([1/6,-1/6,0;1/6,1/3,0;1/6,5/6,0],[1/6,2/3,1/6]),eps);
%!assert(lobatto3b_abc(4),rkinit_abc([1/12,(-1-sqrt(5))/24,(-1+sqrt(5))/24,0;1/12,(25+sqrt(5))/120,(25-13*sqrt(5))/120,0;1/12,(25+13*sqrt(5))/120,(25-sqrt(5))/120,0;1/12,(11-sqrt(5))/24,(11+sqrt(5))/24,0],[1/12,5/12,5/12,1/12]),100*eps);
