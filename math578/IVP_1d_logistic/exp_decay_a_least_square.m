function [rho,C] = exp_decay_a_least_square(a)

abs_a=abs(a);

m = length(abs_a);

Y = log10(abs_a);

k=(0:m-1)';
sum_k=sum(k);
sum_k2=sum(k.^2);

sum_Y=sum(Y);
sum_kY=sum(k.*Y);
denom=m*sum_k2-sum_k^2;
A=(m*sum_kY-sum_k*sum_Y)/denom;
B=(sum_Y*sum_k2-sum_k*sum_kY)/denom;

C = 10^B;
rho = 10^(-A);

%display(['C = ', num2str(C),', rho = ', num2str(rho)])

figure
plot(k,A*k+B,'color',[1 0 0],'linewidth',3)
hold on
plot(Y,'*','color',[0 0 0])

set(gca,'FontSize',20)

axis tight
xlabel('$$k$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$\log_{10}(|a_k|)$$', 'Interpreter', 'latex', 'FontSize', 30)

end






