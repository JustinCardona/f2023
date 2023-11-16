function q4e()
    addpath 'C:\Users\rjust\f2023\math578'
    % Code to produce graphs for A1.4e
    %[T, X, iter] = ODESolve(@(t) 1000 * (t^2 + 1), ...
    %                @(t) 20 * cos(pi * t), ...
    %                @(t) sin(pi * t), ...
    %                2, -1, 0, 1, 2500, "J");
    %[T, X, iter] = ODESolve(@(t) 10^8 * (t^2 + 1), ...
    %                @(t) (1/20) * cos(pi * t), ...
    %                @(t) sin(pi * t), ...
    %                2, -1, 0, 1, 2500, "J");
    [T, X, iter] = ODESolve(@(t) 10^8 * (t^2 + 1), ...
                    @(t) (1/20) * cos(pi * t), ...
                    @(t) sin(pi * t), ...
                    2, -1, 0, 1, 2500, "LU");
    plot(T, X);
    fprintf("%i\n", iter);
end