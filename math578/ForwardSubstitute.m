function y = ForwardSubstitute(L, b)
    % Inputs:
    %   L: A square lower triangular matrix
    %   b: a vector of the length as the columns of U
    % Outputs:
    %   y: A vector that is the solution to Ly = b.
    n = size(L, 1);
    % obtain the first element of y
    y = b / L(1, 1);
    % obtain the other elements of y by iterating downwards in L.
    for k=2:n
        y(k) = (1 / L(k, k)) * (b(k) - L(k, 1:k-1)*y(1:k-1));
    end
end