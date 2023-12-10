addpath '/home/justin/f2023/math578/'
Ns = [20, 40, 80, 160, 320, 640, 1280];
figure
E = zeros(length(Ns), 1);
for k=1:length(Ns)
    x = linspace(-1, 1, Ns(k)).';
    y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
    hold on
    [X, s] = CubicSpline(x, y, [50/676, -50/676]);
    x = linspace(-1, 1, 1+(Ns(k)-1)*10).';
    y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
    E(k) = max(y-s.');
end
x = linspace(-1, 1, 12800).';
y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
plot(x, y)
legend(["True Solution", string(Ns)])
hold off;
fit = fitlm(log(2 ./ (Ns - 1)), log(E)).Coefficients;
plot(log(2 ./ (Ns - 1)), log(E))
fit(2, 'Estimate')

