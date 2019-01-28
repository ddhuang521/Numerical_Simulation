
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

nprint = 20;
iprint = floor(nt/nprint);

% initial solution

umax = .1;

w = umax*exp(-40*(x - L/2).^2);

plot(x,w);
set(gcf,'color','w')
xlabel('x','FontSize',16);
ylabel('u','FontSize',16);
title(['time t = ',num2str(t(1))]);
axis([0 L -2*umax 2*umax]);
pause(0.5);

% Neumann BC

w(1) = (4*w(2)-w(3))/3;
w(nx) = (4*w(nx-1)-w(nx-2))/3;

uexac(1) = 0;
uexac(nx) = 0;

for j=2:nx-1
   
    u(j) = (w(j+1)+w(j-1))/2;
    
end

u(1) = (4*u(2)-u(3))/3;
u(nx) = (4*u(nx-1)-u(nx-2))/3;

% numerical solution

for n=2:nt
    
	for j=2:nx-1

		v(j) = mu^2*u(j+1)+(2-2*mu^2)*u(j)+mu^2*u(j-1)-w(j);
        
		uexac(j) = (umax*exp(-40*(x(j) - L/2 - c*t(n)).^2)+umax*exp(-40*(x(j) - L/2 + c*t(n)).^2))/2;
        
    end
    
    uexac(1) = (4*uexac(2)-uexac(3))/3;
    uexac(nx) = (4*uexac(nx-1)-uexac(nx-2))/3;    

    
    w(1) = (4*w(2)-w(3))/3;
    w(nx) = (4*w(nx-1)-w(nx-2))/3;
        
    u(1) = (4*u(2)-u(3))/3;
    u(nx) = (4*u(nx-1)-u(nx-2))/3;
    
    v(1) = (4*v(2)-v(3))/3;
    v(nx) = (4*v(nx-1)-v(nx-2))/3;
    
	w = u;
	u = v;
	
	if ((n == nt) | (mod(n-1,iprint) == 0))
		plot(x,u,'b-',x,uexac,'r--');
        set(gcf,'color','w')
		xlabel('x','FontSize',16);
		ylabel('u','FontSize',16);
		title(['time t = ',num2str(t(n))]);
        legend('u','uexac')
		axis([0 L -umax umax]);
		pause(0.5);
	end
end





