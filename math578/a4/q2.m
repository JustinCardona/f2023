mu = 2;
L = 0.2;
x0 = [1/100, 1/100];
N = 25;
M = 100;


x = zeros(M, 2);
x(1, :) = x0;

for step=2:M
    a = coeff(x(step - 1), N, L, mu);
    x(step, :) = sum(a, 1);
end
plot(linspace(0, M*L, M), x(:, 1))
hold on
plot(linspace(0, M*L, M), x(:, 2))
legend('x_1', 'x_2')

function a = coeff(x_init, N, L, mu)
    a = zeros(N, 2);
    a(1, :) = x_init;
    for k=2:N
        c = cauchy(cauchy(a(:, 1), a(:, 1)), a(:, 2));
        a(k, 1) = (L / k) * a(k - 1, 2);
        a(k, 2) = (L / k) * mu * (a(k - 1, 2) - c(k-1)) - a(k-1, 1);
    end
end

function c = cauchy(a, b)
    N = size(a);
    c = zeros(N);
    for k=1:N - 1
        c(k) = sum(flip(a(N - k + 1:N)) .* b(1:k));
    end
end