function A = LUDecomposition(A)
    % Inputs:
    %   A: A square matrix
    % Outputs:
    %   A: A square matrix whose upper and lower triangle (below the
    %   diagonal) are the respective U and L such that A = LU.
    n = size(A, 1);
    for k=1:n-1
        for i=k+1:n
            A(i, k) = A(i, k) / A(k, k);
            for j=k+1:n
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
            end
        end
    end
end

