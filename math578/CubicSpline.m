function [X, s] = CubicSpline(x, y, z)
    N = length(x);
    n = N - 1;
    h = (x(N) - x(1)) / n;
    A = diag(4 * ones(1, n+1)) + diag(ones(1,n), 1) + diag(ones(1, n), -1);
    for  i=1:n-1
        u(i+1) = 3/h^2*(y(i+2)-2*y(i+1)+y(i));
    end
    u(1) = 3/h*(y(2) - y(1)) - z(1);
    u(n+1) = z(2) - 3/h*(y(N) - y(N-1));
    u = u';
    % Solving Ac = u from the notes, c 
    w = GaussSeidel(A, u, 1e-12);
    c = [0; w; 0];
    for i = 1:n
        a(i) = y(i);
        d(i) = (c(i+1)-c(i))/(3*h);
        b(i) = (y(i+1)-y(i))/h - h/3*(2*c(i)+c(i+1));
    end
    a(N) = y(N);
    b(1) = z(1);
    b(N) = z(2);
    m = 10;
    h = h/m;
    X=x(1):h:x(N);
    for i = 1:n
        for j = m*(i-1)+1:m*i
            s(j) = a(i) + b(i)*(X(j)-x(i)) + c(i)*(X(j)-x(i))^2 + d(i)* (X(j)-x(i))^3;
        end
    end
    s(m*n+1) = y(N);
    plot(X,s)
end

