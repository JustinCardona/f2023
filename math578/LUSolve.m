function x = LUSolve(A, b)
    % Inputs:
    %   A: A square matrix
    %   b: a vector of the length as the columns of A
    % Outputs:
    %   x: A vector that is the solution to Ax = b done by LU Decomposition.
    A = LUDecomposition(A);
    % Perform the LU Decomposition of A
    L = tril(A, -1) + eye(size(A, 1));
    U = triu(A);
    % Solve Ly = b and Ux = y in order to solve Ax = b
    x = BackwardSubstitute(U, ForwardSubstitute(L, b));
end

