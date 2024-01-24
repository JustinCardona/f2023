function [] = plot_a_taylor(a,x0,nu)

N=length(a)-1;

t=(-nu:.001:nu);

m=length(t);
u=zeros(1,m);

for k=1:m
    u(k)=sum(a.*t(k).^(0:N)');
end
    

%color = 'red';
color = 'blue';
%color = 'green';
plot(t,u,'Linewidth',3,'color',color)
hold on

scatter(0,x0,100,'black','filled')

plot(0,x0,'*','Linewidth',3,'color',[0 0 0])

axis tight

set(gca,'FontSize',20)

xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$x(t)$$', 'Interpreter', 'latex', 'FontSize', 30)

end

    