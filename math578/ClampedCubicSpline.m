function [a, b, c, d] = ClampedCubicSpline(x, y, z)
    n = size(x, 1);
    M = zeros(n, n);
    h = x(2:n) - x(1:n-1);
    f = y(2:n) - y(1:n-1) ./ h;
    for k=1:n-1
        M(k, k+1) = h(k);
        M(k+1, k) = h(k);
    end
    for k=2:n-1
        M(k, k) = 2*(h(k) + h(k-1));
    end
    M(1, 1) = 2*h(1);
    M(n, n) = 2*h(n-1);
    m = zeros(n, 1);
    m(1) = f(1) - z(1);
    m(2:n-1) = 3 * (f(2:n-1) - f(1:n-2));    
    m(n) = 3 * (z(2) - f(n-1));
    c = M\m;%GaussSeidel(M, m, 1e-14);
    d = (c(2:n) - c(1:n-1)) ./ (3*h);
    d(n) = (c(n) - c(n-1)) / (3*h(n-1));
    b = f - (h/3) .* (c(2:n) + 2*c(1:n-1));
    a = y;
    b(1) = z(1);
    b(n) = z(2);
    d
    N = 10;
    y_new = zeros(N*n, 1);
    for k=1:n-1
        domain = linspace(x(k), x(k+1), N).';
        y_new(N*(k-1)+1:N*k) = arrayfun(@(t) a(k) + b(k)*(t-x(k)) + c(k)*(t-x(k))^2 + d(k)*(t-x(k))^3, domain);
    end
    plot(linspace(min(x), max(x), N*n).', y_new.')
end

