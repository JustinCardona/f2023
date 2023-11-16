function c = LeastSquaresQuadratic(x, y)
    % Inputs:
    %   x: Vector of independant of the data
    %   y: Vector of dependant of the data
    % Outputs:
    %   c: Vector representing the coefficients of the best fitting
    %   quadratic to the data
    % Form M from x
    M = [ones(size(x)) x x.*x];

    % Solve M^TMc = M^Ty
    c = LUSolve(M.' * M, M.' * y);
end

