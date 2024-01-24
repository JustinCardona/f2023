addpath '/home/justin/f2023/math578/intlab/intlab_V11'

% Choose initial condition
x0 = 4/5;

% Choose the Taylor order of the finite dimensional reduction
N = 50;

% Compute Taylor coefficients
a = zeros(N+1,1);
a(1) = x0;
a(2) = -8/25;
a(3) = -64/125;

for k = 4:N+1
    a(k) = (1/k)*(a(k-3) - 0.5*sum(a(1:k-1).*a(k-1:-1:1)));
    %a(k) = (1/(k-1))*(-a(k-1)+sum(a(1:k-1).*a(k-1:-1:1)));
end

[rho,C] = exp_decay_a_least_square(a);
disp(['Geometric decay = ',num2str(rho)])

% Choose the size nu of the domain (-nu,nu) on which we attempt the proof
%nu = floor(rho);
nu = 1;
[I] = int_rad_poly(x0,nu,N);
disp(['I = [',num2str((I(1))),' , ',num2str((I(2))),']'])

figure
plot_a_taylor(a,x0,nu)