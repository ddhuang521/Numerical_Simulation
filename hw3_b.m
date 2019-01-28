
% Solve the 1D wave equation by the upwind scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; format long;

L = 10;
T = 20;

c = 1;
dx = 0.01;
dt = 0.01;

mu = c*dt/dx;

x = 0:dx:L;
t = 0:dt:T;

nx = length(x);
nt = length(t);

% display solution

nprint = 5;
iprint = floor(nt/nprint);

% initial solution

umax = 0.1;

w = umax*exp(-40*(x - L/2).^2);

plot(x,w);
set(gcf,'color','w')
xlabel('x','FontSize',16);
ylabel('u','FontSize',16);
title(['time t = ',num2str(t(1))]);
axis([0 L -umax umax]);
pause(0.5);

% Dirichlet BC

w(1) = 0;
w(nx) = 0;
u(1) = 0;
u(nx) = 0;
v(1) = u(1);
v(nx) = u(nx);
uexac(1) = u(1);
uexac(nx) = u(nx);

for j = 2:nx-1
   
    u(j) = (w(j+1)+w(j-1))/2;
    
end

% numerical solution

for n=2:nt
    
	for j=2:nx-1
        
		v(j) = mu^2*u(j+1)+(2-2*mu^2)*u(j)+mu^2*u(j-1)-w(j);
        
		uexac(j) = (umax*exp(-40*(x(j) - L/2 - c*t(n)).^2)+umax*exp(-40*(x(j) - L/2 + c*t(n)).^2))/2;
        
    end
    
	w = u;
	u = v;
	
	if ((n == nt) | (mod(n-1,iprint) == 0))
		plot(x,u,'-',x,uexac,'--');
        set(gcf,'color','w')
		xlabel('x','FontSize',16);
		ylabel('u','FontSize',16);
		title(['time t = ',num2str(t(n))]);
		axis([0 L -umax umax]);
        legend('u','uexac')
		pause(0.5);
	end
end





