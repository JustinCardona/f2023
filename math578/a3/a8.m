eps = 0.03;
init = 0.1;
t_max = 20;
addpath '/home/justin/f2023/math578/'
figure; hold on
[t, y] = Euler(@(x)f(x), 0, t_max, init, 250);
plot(t, y);
[t, y] = TrapezoidEuler(@(x)f(x), 0, t_max, init, 500);
plot(t, y);

function Y = f(y)
    Y = y^2 - eps * y^3
end
