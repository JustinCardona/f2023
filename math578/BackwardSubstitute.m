function x = BackwardSubstitute(U, b)
    % Inputs:
    %   U: A square upper triangular matrix
    %   b: a vector of the length as the columns of U
    % Outputs:
    %   x: A vector that is the solution to Ux = b.
    n = size(U, 1);
    % obtain the last element of x
    x = b / U(n, n);
    % obtain the other elements of x by iterating upwards in U.
    for k=n-1:-1:1
        x(k) = (1 / U(k, k)) * (b(k) - U(k, k+1:n) * x(k+1:n));
    end
end