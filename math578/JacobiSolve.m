function [x, iter] = JacobiSolve(A, b)
    % Inputs:
    %   A: A square matrix
    %   b: a vector of the length as the columns of A
    % Outputs:
    %   x: A vector that is the solution to Ax = b done by the Jacobi method.
    n = size(b, 1);
    x = rand(n, 1);
    iter = 0;
    D_inverse = diag(1 ./ diag(A));
    B_J = D_inverse * (-tril(A, -1) - triu(A, 1));
    B_J = D_inverse * (-tril(A, -1) - triu(A, 1));
    while norm(A*x - b, inf) > 1e-6
        x = B_J*x + D_inverse*b;
        iter = iter + 1;
    end
end

