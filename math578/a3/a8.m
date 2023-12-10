epsilon = 0.03;
init = 0.1;
t_max = 20;
addpath '/home/justin/f2023/math578/'
figure
hold on
[t, y] = Euler(@(x)f(x, epsilon), 0, t_max, init, 250);
plot(t, y);
hold on
[t, y] = TrapezoidEuler(@(x)f(x, epsilon), 0, t_max, init, 250);
plot(t, y);
legend('Euler', 'Trapezoid')

function Y = f(y, epsilon)
    Y = y^2*(1 - epsilon * y);
end

