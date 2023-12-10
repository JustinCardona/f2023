function x = GaussSeidel(A, b, tol)
    E = -tril(A, -1);
    F = -triu(A, 1);
    D = diag(diag(A));
    B = (D-E)\F;
    add = (D-E)\b;
    x = zeros(size(b));
    while norm(A*x-b, 2)>tol
        x = B*x + add;
    end
end