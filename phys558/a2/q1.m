syms a b c A B C;
S = sqrt(sin(A)^2*sin(B)^2 - (cos(C) - cos(A)*cos(B))^2);
V = simplify([a 0 0; b*cos(A) b*sin(A) 0; c*cos(B) c*(cos(C)-cos(A)*cos(B))/sin(A) c/(sin(A))*S]);
U = simplify(2 * pi * inv(simplify(V)));
W = simplify(U.' * U);
latex(W(1:3))
latex(W(4:6))
latex(W(7:9))