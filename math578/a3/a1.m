syms f(x1, x2, x3, x4, lambda) Delta(lambda) x1 x2 x3 x4 lambda
Delta(lambda) = 9.32 * log(11/10 - lambda/500);
f(x1, x2, x3, x4, lambda) = [
    -20 - 120*x2^2*x3*(x1 - 25.1*log((23/1350) * (550 - lambda))) - 36*x4^4*(x1 - 25.1*log(lambda/400)) - 0.3*(x1 + 24.3) ...
    (9/25)*(1-x2)*(x1-Delta(lambda)+35)/(1-exp(-(x1-Delta(lambda)+35)/10)) - (72/5)*x2*exp(-(x1-Delta(lambda)+60)/18) ...
    (63/250)*(1-x3)*exp(-(x1-Delta(lambda)+60)/20) - (18/5)*x3/(exp(-(x1-Delta(lambda)+30)/10)+1) ...
    (9/250)*(1-x4)*(x1-Delta(lambda)+50)/(1 - (1+exp(-(x1-Delta(lambda)+50)/10))) - (9/20)*x4*exp(-(x1-Delta(lambda)+60)/80)...
];

X0 = [-90.9665063496862; 0.000846010922010574; 0.993896751841571; 0.0341288422401492; 11];    
Ds = 0.05;
steps = 25000;
h = Ds / 1e3;
cont = zeros(steps, 5);
DXf = [
        (f(X0(1)+h/2, X0(2), X0(3), X0(4), X0(5)) - f(X0(1)-h/2, X0(2), X0(3), X0(4), X0(5)))/h;...
        (f(X0(1), X0(2)+h/2, X0(3), X0(4), X0(5)) - f(X0(1), X0(2)-h/2, X0(3), X0(4), X0(5)))/h;...
        (f(X0(1), X0(2), X0(3)+h/2, X0(4), X0(5)) - f(X0(1), X0(2), X0(3)-h/2, X0(4), X0(5)))/h;...
        (f(X0(1), X0(2), X0(3), X0(4)+h/2, X0(5)) - f(X0(1), X0(2), X0(3), X0(4)-h/2, X0(5)))/h;...
        (f(X0(1), X0(2), X0(3), X0(4), X0(5)+h/2) - f(X0(1), X0(2), X0(3), X0(4), X0(5)-h/2))/h;...
    ];
DXf = double(DXf).';
for k=1:steps
    [U, S, V] = svd(DXf);
    X0dotNew = V(:,5) / norm(V(:,5), 2);
    if dot(X0dotNew, X0dot) < 0
        X0dot = -X0dotNew;
    else
        X0dot = X0dotNew;
    end

    X1hat = X0 + Ds*X0dot;
    Dinvf = inv(cat(1, X0dot.', DXf));
    
    %Newton Method
    X0 = X1hat;
    e = cat(1, dot((X0 - X1hat), X0dot), double(f(X0(1), X0(2), X0(3), X0(4), X0(5))).');
    tol = 1e-12;
    while norm(e) > tol
        DXf = [
        (f(X0(1)+h/2, X0(2), X0(3), X0(4), X0(5)) - f(X0(1)-h/2, X0(2), X0(3), X0(4), X0(5)))/h;...
        (f(X0(1), X0(2)+h/2, X0(3), X0(4), X0(5)) - f(X0(1), X0(2)-h/2, X0(3), X0(4), X0(5)))/h;...
        (f(X0(1), X0(2), X0(3)+h/2, X0(4), X0(5)) - f(X0(1), X0(2), X0(3)-h/2, X0(4), X0(5)))/h;...
        (f(X0(1), X0(2), X0(3), X0(4)+h/2, X0(5)) - f(X0(1), X0(2), X0(3), X0(4)-h/2, X0(5)))/h;...
        (f(X0(1), X0(2), X0(3), X0(4), X0(5)+h/2) - f(X0(1), X0(2), X0(3), X0(4), X0(5)-h/2))/h;...
        ];
        DXf = double(DXf).';
        Dinvf = inv(cat(1, X0dot.', DXf));
        X0 = X0 - Dinvf*e;
        e = cat(1, dot((X0 - X1hat), X0dot), double(f(X0(1), X0(2), X0(3), X0(4), X0(5))).');
    end
    cont(k,:) = X0;
end

plot(cont(:, 1), cont(:, 5));
