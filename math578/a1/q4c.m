function q4c()
    % Code to produce graphs for A1.4
    n = arrayfun(@(x) 10 * 2^x, 0:5);    
    for k=0:length(n)-1
        [T, X, ~] = ODESolve(@(t) 100 * (t^2 + 1), ...
                        @(t) 20 * cos(pi * t), ...
                        @(t) sin(pi * t), ...
                        2, -1, 0, 1, n(k+1), "LU");
        plot(T, X, 'DisplayName',"n = " + k);
        hold on;
    end
    legend()
end