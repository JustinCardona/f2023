function [t,x] = ImprovedEuler(f,a,b,x0,n, d)
    h = (b-a)/n;
    t = linspace(a,b,n+1);
    x = zeros(n, d);
    x(1, :) = x0;
    for i = 1:n
        x_eul = x(i, :) + h*f(x(i, :));
        x(i+1, :) = x(i, :) + h*(f(x(i, :))+f(x_eul))/2;
    end