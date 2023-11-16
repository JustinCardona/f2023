function cont = PseudoArcLengthConituation(X0, f, D, steps)
    cont = zeros(steps, 1);
    for k=1:steps        
        DXf(x1, x2, x3, x4, lambda) = [diff(f, x1); diff(f, x2); diff(f, x3); diff(f, x4); diff(f, lambda)];
        X0bar = double((DXf(X0(1), X0(2), X0(3), X0(4), X0(5)))).' * X0
    end

end

