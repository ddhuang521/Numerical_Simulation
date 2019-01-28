
% Solve the 1D wave equation by the upwind scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wave_upwind

clc; clear all; format long;

L = 10;
T = 20;

c = input('Enter wave speed (> 0) ==> ');
dx = input('Enter grid size dx ==> ');

dt = input('Enter time step dt ==> ');

mu = c*dt/dx;

x = 0:dx:L;
t = 0:dt:T;

nx = length(x);
nt = length(t);

% display solution

nprint = 50;
iprint = floor(nt/nprint);

% initial solution

umax = 0.1;

u = umax*exp(-40*(x - L/4).^2);

plot(x,u);
xlabel('x','FontSize',16);
ylabel('u','FontSize',16);
title(['time t = ',num2str(t(1))]);
axis([0 L -2*umax 2*umax]);
pause(0.5);

% boundary conditions

u(1) = 0;
u(nx) = 0;
v(1) = u(1);
v(nx) = u(nx);
uexac(1) = u(1);
uexac(nx) = u(nx);

for n=2:nt
	
	for j=2:nx-1
		v(j) = (1 - mu)*u(j) + mu*u(j-1);
		uexac(j) = umax*exp(-40*(x(j) - L/4 - c*t(n)).^2);
	end
	
	u = v;
	
	if ((n == nt) | (mod(n,iprint) == 0))
		plot(x,u,'-',x,uexac,'--');
		xlabel('x','FontSize',16);
		ylabel('u','FontSize',16);
		title(['time t = ',num2str(t(n))]);
		axis([0 L -2*umax 2*umax]);
		pause(0.5);
	end
end





