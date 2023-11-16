addpath '/home/justin/f2023/math578/'
Ns = [20, 40, 80, 160, 320, 640, 1280];
figure;hold on
for k=1:size(Ns)
    x = linspace(-1, 1, Ns(k)).';
    y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
    [a, b, c, d] = CubicSpline(x, y, [50/676, -50/676]);
end
x = linspace(-1, 1, 2560).';
y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
plot(x, y)