addpath '/home/justin/f2023/math578/'

[t, x] = Euler(@(x) f(x), 0, 30, [pi/4; 0], 1000, 2);
[ti, xi] = ImprovedEuler(@(x) f(x), 0, 30, [pi/4; 0], 1000, 2);
Ei = (xi(:, 2).^2) / 2 - cos(xi(:, 1));
E = (x(:, 2).^2) / 2 - cos(x(:, 1));


%plot(t, x(:, 1))
%hold on
%plot(t, xi(:, 1))
%legend('Euler \theta', 'Improved Euler \theta')

%plot(x(:, 1), x(:, 2))
%hold on
%plot(x(:, 1), xi(:, 2))
%legend('Euler', 'Improved Euler')


plot(t, E)
hold on
plot(t, Ei)
legend('Euler Energy', 'Improved Euler Energy')


function d = f(x)
    d = [x(2), -sin(x(1))];
end