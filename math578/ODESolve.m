function [time, x, iter] = ODESolve(Q_fn, P_fn, B_fn, alpha, beta, t, T, n, method)
    % Inputs:
    %   Q_fn: coefficient function of the first derivative
    %   P_fn: coefficient function of the function
    %   B_fn: added function
    %   a: initial value scalar
    %   b: final value scalar
    %   t: initial time
    %   T: final time
    %   n: number of finite elements
    %   m: method of Solving to be used "LU" for LU and "J" for Jacobi
    % Outputs:
    %   time: The vector representation of the time
    %   x: The vector representation of the time series solution to the ODE.
    %   iter: the number of iteration of the algorithm used
    time = linspace(t, T, n);
    Q = arrayfun(Q_fn, time);
    P = arrayfun(P_fn, time);
    B = arrayfun(B_fn, time);
    A = zeros(n);
    h = (T - t) / n;
    for k=2:n-1
        A(k-1, k) = -h * P(k) / 2 - 1;
        A(k, k) = h^2 * Q(k) + 2;
        A(k+1, k) = h * P(k) / 2 - 1;
    end
    A(1, 1) = h^2 * Q(1) + 2;
    A(2, 1) = h * P(1) / 2 - 1;
    A(n, n) = h^2 * Q(n) + 2;
    A(n-1, n) = -h * P(n) / 2 - 1;
    b = -h^2 * B;
    b(1) = b(1) + (1 + h * P(1) / 2) * alpha;
    b(n) = b(n) + (1 - h * P(n) / 2) * beta;
    switch method 
        case "LU"
            x = LUSolve(A, b.');
            iter = 1;
        otherwise
            [x, iter] = JacobiSolve(A, b.');
    end
end