clc;clear all;

% part A

Np=[2:2:100]; % grid points pool

for a=1:length(Np)
    
    N=Np(a);

    n=[0:1:N-1];
    dx=2*pi/N;
    dx1(a)=dx;
    x=-pi+dx.*n;

    u=exp(sin(x)); 
    u1=exp(sin(x)).*cos(x);
    u2=exp(sin(x)).*cos(x).^2-sin(x).*exp(sin(x));

    k=[-N/2:1:N/2-1];

    for j=1:N
    
    uk(j)=sum(u.*exp(-1i.*k(j).*x))/N;
    
    end

    % part A 1st order derivative
    
    % DFT
    
    for m=1:N

    uN1(m)=sum(1i.*k.*uk.*exp(1i.*k.*x(m)));

    end

    err1(a)=max(abs(uN1-u1));
    
    % finite difference
    
    ud1(1)=(u(2)-u(N))/(2*dx);
    ud1(N)=(u(1)-u(N-1))/(2*dx);
    
    for j=2:N-1
        
    ud1(j)=(u(j+1)-u(j-1))/(2*dx);
    
    end
    
    
    err2(a)=max(abs(ud1-u1));
    
    % part B 2nd order derivative
    
    % fft
    
    for m=1:N

    uN2(m)=sum(-k.^2.*uk.*exp(1i.*k.*x(m)));

    end
    
    err3(a)=max(abs(uN2-u2));
    
    % finite difference    
    
    ud2(1)=(u(2)+u(N)-2*u(1))/(dx^2);
    ud2(N)=(u(1)+u(N-1)-2*u(N))/(dx^2);
    
    for j=2:N-1
        
    ud2(j)=(u(j+1)+u(j-1)-2*u(j))/(dx^2);
    
    end
    
    
    err4(a)=max(abs(ud2-u2));
    

end

figure;
loglog(dx1,err1,dx1,err2,'ro');
set(gcf,'color','w');
title('1st order derivative logerr vs logdx','Fontsize',16);
xlabel('log dx','Fontsize',16);
ylabel('log err','Fontsize',16);
legend('DFT','FD');

figure;
loglog(dx1,err3,dx1,err4,'ro');
set(gcf,'color','w');
title('2nd order derivative logerr vs logdx','Fontsize',16);
xlabel('log dx','Fontsize',16);
ylabel('log err','Fontsize',16);
legend('DFT','FD');

% figure;
% plot(x,u1);
% hold on;
% plot(x,uN1,'ro');
% hold on;
% plot(x,ud1,'c*');
% set(gcf,'color','w');
% title('1st order derivative of u','Fontsize',16);
% xlabel('x','Fontsize',16);
% ylabel('u1','Fontsize',16);
% legend('analytic','DFT','FD');
% 
% figure;
% plot(x,u2);
% hold on;
% plot(x,uN2,'ro');
% hold on;
% plot(x,ud2,'c*');
% set(gcf,'color','w');
% title('2nd order derivative of u','Fontsize',16);
% xlabel('x','Fontsize',16);
% ylabel('u2','Fontsize',16);
% legend('analytic','DFT','FD');



% figure;
% semilogy(Np,err1);
% hold on;
% semilogy(Np,err2,'ro');
% set(gcf,'color','w');
% title('1st order derivative errors comparison','Fontsize',16);
% xlabel('N(grid points)','Fontsize',16);
% ylabel('log errors','Fontsize',16);
% legend('DFT','FD','Fontsize',16);
% 
% figure;
% semilogy(Np,err3);
% hold on;
% semilogy(Np,err4,'ro');
% set(gcf,'color','w');
% title('2nd order derivative errors','Fontsize',16);
% xlabel('N(grid points)','Fontsize',16);
% ylabel('log errors','Fontsize',16);
% legend('DFT','FD','Fontsize',16);



