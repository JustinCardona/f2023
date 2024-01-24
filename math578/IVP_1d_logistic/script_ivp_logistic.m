addpath '/home/justin/f2023/math578/IVP_1d_logistic/'
clear
clc
close all

% Choose initial condition
x0 = 5/6;

% Choose the Taylor order of the finite dimensional reduction
N = 200;

% Compute Taylor coefficients
a = zeros(N+1,1);
a(1) = x0;

for k = 2:N+1
    a(k) = (1/(k-1))*(-a(k-1)+sum(a(1:k-1).*a(k-1:-1:1)));
end

[rho,C] = exp_decay_a_least_square(a);
disp(['Geometric decay = ',num2str(rho)])

% Choose the size nu of the domain (-nu,nu) on which we attempt the proof
%nu = floor(rho);
nu = 2.9;
[I] = int_rad_poly(x0,nu,N);
disp(['I = [',num2str((I(1))),' , ',num2str((I(2))),']'])

figure
plot_a_taylor(a,x0,nu)

