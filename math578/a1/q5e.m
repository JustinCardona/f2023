function c = q5e()
    % code to produce graphs and error margins for 1.5e
    % form points
    addpath 'C:\Users\rjust\f2023\math578'
    x = [5, 5.5, 6.5, 8, 8.5, 10.8, 11.5, 13.7, 14.5, 15.9];
    y = [1, 4, 7, 8, 9.5, 9.2, 9, 6, 3, 1];

    % fitting
    c = LeastSquaresQuadratic(x.', y.');    
    fit_fn = @(X) c(1) + c(2)*X + c(3)*X^2;

    % calculate residuals
    fprintf("%d\n", sqrt(sum((y - arrayfun(fit_fn, x)).^2)));

    % plot
    scatter(x, y)
    hold on;
    t = linspace(min(x), max(x));
    plot(t, arrayfun(fit_fn, t));
end

