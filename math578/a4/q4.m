A = [1, 2, 3;
     2, 3, 4;
     3, 4, 5];
[Q, R] = QRDecompose(A)

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