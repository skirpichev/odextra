function abc = collocation_abc(P)
abc.s = length(P) - 1;
abc.c = sort(roots(P));
V = vander(abc.c);
V = V(:,abc.s:-1:1);
J = diag(1 ./ (1:abc.s)) / V;
abc.a = diag(abc.c) * V * J;
abc.b = ones(1,abc.s) * J;
