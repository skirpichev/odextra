% Shifted Nth-order Legendre polynomials normalized in the [0,1].

function P = legendre_coeff(s)
Pm1 = [];
P = 1;
for m = 1:s
  A(1) = m;
  A(2) = -2 * A(1) + 1;
  A(3) = -2 * A(2);
  A(4) = -A(1) + 1;
  PP1 = A(2) * [0,P] + A(3) * [P,0] + A(4) * [0,0,Pm1];
  PP1 = PP1 / A(1);
  Pm1 = P;
  P = PP1;
end

%!assert(legendre_coeff(0),1);
%!assert(legendre_coeff(1),[2,-1]);
%!assert(legendre_coeff(2),[6,-6,1]);
%!assert(legendre_coeff(7),[3432,-12012,16632,-11550,4200,-756,56,-1]);
