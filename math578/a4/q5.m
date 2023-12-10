addpath '/home/justin/f2023/math578/'

A = hilb(9);
%[lambda_1, u_1, count_1] = DPM(A, 1e-14)
[lambda_9, u_9, count_9] = IPM(A, 1e-14)
%[lambda, count] = QRIM(A, 1e-14)


function [v, q, count] = DPM(A, tol)
    q0 = rand(9, 1, "like", 1i);
    q = q0 / norm(q0, 2);

    z = A*q;
    q = z / norm(z, 2);
    v = conj(q).' * A * q;
    count = 1;
    while norm(z - v*q, 2) > tol
        z = A*q;
        q = z / norm(z, 2);
        v = conj(q).' * A * q;
        count = count + 1;
    end
end


function [v, q, count] = IPM(A, tol)
    q0 = rand(9, 1, "like", 1i);
    q = q0 / norm(q0, 2);

    M = LUDecomposition(A);
    L = tril(M, -1) + eye(size(M, 1));
    U = triu(M);

    z = BackwardSubstitute(U, ForwardSubstitute(L, q));
    q = z / norm(z, 2);
    v = conj(q).' * A * q;
    err = norm(z - v*q, 2)
    count = 1;
    while err > tol
        z = BackwardSubstitute(U, ForwardSubstitute(L, q))
        q = z / norm(z, 2);
        v = conj(q).' * A * q;
        count = count + 1;
        err = norm(z - v*q, 2);
    end
end


function [lambda, count] = QRIM(A, tol)
    [Q, R] = QRDecompose(A);
    T = R * Q
    r = norm(diag(T, 1), 2)
    count = 0;
    while r > tol
        [Q, R] = QRDecompose(T);
        T = R * Q;
        r = norm(diag(T, 1), 2)
        count = count + 1;
    end
    lambda = diag(T);
end

function [Q, R] = QRDecompose(A)
    [~,n] = size(A);
    Q = A;
    R = zeros(n);    
    for j = 1:n
        R(j,j) = norm(Q(:,j));
        Q(:,j) = Q(:,j)/R(j,j);
        R(j,j+1:n) = Q(:,j)'*Q(:,j+1:n);
        Q(:,j+1:n) = Q(:,j+1:n) - Q(:,j)*R(j,j+1:n);
    end
end