clc;clear all;

Np=[2:2:100];

for w=1:length(Np)
    
    N=Np(w);
    
n=[0:1:N];
dx=pi/N;
x=cos(dx.*n);

% analytic solutions

u=1./(1+x.^2);
u1=-2.*x./((1+x.^2).^2); % u'
u2=8*x.^2./((1+x.^2).^3)-2./((1+x.^2).^2); % u''

% matrix DN construction

D=zeros(N+1,N+1);

c=zeros(1,N+1);
c(1)=2;
c(N+1)=2;
c(2:N)=1;


for l=1:N+1
   for m=1:N+1 
    
       if l~=m
       
           D(l,m)=c(l)*(-1)^(l+m)/(c(m)*(x(l)-x(m)));
          
       else
           
       end
   end
end


for l=1:N+1
    
    D(l,l)=-sum(D(l,:)); 
    
end

uN1=(D*u')'; % 1st order derivative
uN2=(D*D*u')'; % 2nd order derivative


err1(w)=max(abs(uN1-u1)); % max norm
err2(w)=max(abs(uN2-u2));

end

figure;
semilogy(Np,err1);
set(gcf,'color','w');
title('1st order derivative error','FontSize',16);
xlabel('N(grid points)','FontSize',16);
ylabel('log error','FontSize',16);

figure;
semilogy(Np,err2);
set(gcf,'color','w');
title('2nd order derivative error','FontSize',16);
xlabel('N(grid points)','FontSize',16);
ylabel('log error','FontSize',16);








