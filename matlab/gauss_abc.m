function abc = gauss_abc(s)
abc = collocation_abc(legendre_coeff(s));

%!assert(gauss_abc(1),rkinit_abc(0.5,1));
%!assert(gauss_abc(2),rkinit_abc((1 + [0,-2;2,0]/sqrt(3))/4,[1,1]/2),eps);
