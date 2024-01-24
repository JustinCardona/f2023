function I = int_rad_poly(x0,nu,N)

inu = intval(nu);
iN = intval(N);

a = intval(zeros(N+1,1));
a(1) = intval(x0);

for k = 2:N+1
    a(k) = (1/(k-1))*(-a(k-1)+sum(a(1:k-1).*a(k-1:-1:1)));
end

norm_a = sum(inu.^(0:N)'.*abs(a));

%%%%%%%%%%%
%%%% Y0 %%%%
%%%%%%%%%%%

a = [a;zeros(N,1)];
s0 = intval(0);
for k = N:2*N
    s0 = s0 + (1/(k+1))*abs(-a(k+1)+sum(a(1:k+1).*a(k+1:-1:1)))*inu^(k+1);
end

Y0 = sup(s0);

%%%%%%%%%%%%%%
%%%% Z(r) %%%%
%%%%%%%%%%%%%%

Z1 = sup((1+2*norm_a)*inu/(iN+1));
Z2 = sup(2*inu/(iN+1));

fprintf('\n')
display(['Z2 = ',num2str(Z2),', Z1 = ',num2str(Z1),', Y0 = ',num2str(Y0)])
fprintf('\n')

I = sort(roots([Z2 Z1-1 Y0]))';

end
