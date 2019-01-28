clc; clear all; format long;

L = 50;
T = 40;

c = 1;
dx = 0.01;
dt = 0.012;

mu = c*dt/dx;

x = 0:dx:L;
t = 0:dt:T;

nx = length(x);
nt = length(t);

% display solution

nprint = 3;
iprint = floor(nt/nprint);

% Dirichlet BC

u(1:nx) = 0;
u((nx-1)*3/50+1:(nx-1)*5/50+1)=0.1;


% numerical solution

for n=2:nt
    
    
	for j=2:nx-1
        
		v(j) = ((1-mu)*u(j+1)+(1+mu)*u(j-1))/2;
        
    end
    
    v(1)=0;
    v(nx)=0;
    u = v;
	
	if ((n == nt) | (mod(n-1,iprint) == 0))
		plot(x,u,'-');
        hold on;
        set(gcf,'color','w')
		xlabel('x','FontSize',16);
		ylabel('u','FontSize',16);
		title('dt=0.012');
		axis([0 L -0.5 0.5]);
		pause(0.5);
	end
end





