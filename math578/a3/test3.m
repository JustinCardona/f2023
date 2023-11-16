clear; clc;
x = linspace(-1, 1, 2560).';
y = arrayfun(@(t) 1 / (1 + 25*t^2), x);
N = length(x); %number of points 
n =N -1; % nmber of subintervals 
h = (x(N)-x(1))/n;
Trid = diag(4*ones(1,n-1))+ diag(ones(1,n-2),1)+diag(ones(1,n-2),-1);
for  i = 1:n-1
    z(i) = 6/h^2*(y(i+2)-2*y(i+1)+y(i));
end
z = z';
w = inv(Trid)*z;
sigma = [0; w; 0];
for i = 1:n
    d(i) = y(i);
    b(i) = sigma(i)/2;
    a(i) = (sigma(i+1)-sigma(i))/(6*h);
    c(i) = (y(i+1)-y(i))/h-h/6*(2*sigma(i)+sigma(i+1));
end
r = 4; %number of subsubintervals
hh = h/r; %step size of subsubintervals
m=x(1): hh:x(N);
for i =1:n
    for j = r* (i-1)+1:r*i
        s(j) = a(i)* (m(j)-x(i))^3+b(i)*(m(j)-x(i))^2+c(i)*(m(j)-x(i))+ d(i);
    end
end
s(r*n+1) =y(N);
plot(x, y, 'o')
hold on
plot(m,s,'-x')
hold off