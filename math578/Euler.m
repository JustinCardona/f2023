function [t,x] = Euler(f,a,b,x0,n)
h = (b-a)/n;
t = linspace(a,b,n+1);
x(1) = x0;
for i = 1:n
    x(i+1) = x(i) + h*f(x(i));
end


