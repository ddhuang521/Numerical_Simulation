
% Solve the 1D heat equation by explicit Euler method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; format long;

L = 1;
T = 0.5;

dx = input('Enter grid size dx ==> ');
dt = input('Enter time step dt ==> ');
the= input('Enter the theta parameter ==>');

% dt=dx^2/2;

mu = dt/(dx^2);

x = 0:dx:L;
t = 0:dt:T;

nx = length(x);
nt = length(t);

% display solution

nprint = 50;
iprint = floor(nt/nprint);

% boundary conditions

f0 = exp(-(pi^2)*t);
f1 = -exp(-(pi^2)*t);

% initial solution

u = cos(pi*x);

plot(x,u,'-',x,u,'--');
xlabel('x','FontSize',16);
ylabel('u','FontSize',16);
title(['time t = ',num2str(t(1))]);
legend('numerical','exact');
axis([0 1 -1 1]);
pause(0.5);

% time integration

err = 0;    % maximum error relative to the exact solution
ter = 0;

trun = cputime;

A(1,1)=1;
A(nx,nx)=1;
A(2:nx-1,2:nx-1) = diag((1+2*mu*the)*ones(1,nx-2),0) + diag((-mu*the)*ones(1,nx-3),1) + diag((-mu*the)*ones(1,nx-3),-1);

for n=2:nt-1
    
    u(1) = f0(n);
	u(nx) = f1(n);
    b(1) = f0(n+1);
    b(nx) = f1(n+1);
	uexac(1) = u(1);
	uexac(nx) = u(nx);
    
    
	for j=2:nx-1
        
		b(j) = (1-the)*mu*u(j+1) + (1-2*mu*(1-the))*u(j) - (1-the)*mu*u(j-1);
        
		uexac(j) = cos(pi*x(j))*exp(-(pi^2)*t(n));    % exact solution
        
	end
	
	 u = sparse(A)\b';
    
   
	if ((n == nt) | (mod(n,iprint) == 0))
	    err = [err max(abs(u' - uexac))];
		ter = [ter t(n)];
		plot(x,u,'-',x,uexac,'--');
		xlabel('x','FontSize',16);
		ylabel('u','FontSize',16);
		title(['time t = ',num2str(t(n))]);
		legend('numerical','exact');
		axis([0 1 -1 1]);
		pause(0.5);
	end
end

trun = cputime - trun

figure;

plot(ter,err,'o-');
xlabel('t','FontSize',16);
ylabel('max error','FontSize',16);


