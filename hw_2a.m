clc; clear all; format long;

L = 1;
N = 51;
M = N*N;
J = N+1;

r = linspace(0,L,N);   % x-coordinate
t = linspace(0,pi/2,N);                % y-coordinate

dr = r(2) - r(1);     % grid spacing
dr2 = dr*dr;
dt = t(2) - t(1);
dt2 = dt*dt;

% label the grid points: boundary or interior

for j=1:N
	for i=1:N
		k = i + N*(j - 1);                    % (i,j)    
		ityp(k) = 0;                          % interior point
		if ((i == 1) & (j >= 1) & (j < N))    % r=0
		    ityp(k) = 1;
		end
		if ((i == N) & (j > 1) & (j <= N))    % r=1
		    ityp(k) = 2;
		end
		if ((i > 1) & (i <= N) & (j == 1))    % theta=0
		    ityp(k) = 3;
		end
		if ((i >= 1) & (i < N) & (j == N))    % theta=pi/2
		    ityp(k) = 4;
		end
	end
end

% boundary conditions

% fa = zeros(N,1);
fb = sin(2.*t).^2./4;
%fc = sin(pi*x2);
%fd = exp(pi)*sin(pi*x2) + 0.5*x2.^2;

% set up the linear system

A = zeros(M,M);
b = zeros(M,1);

for j=1:N
    
	for i=1:N
        
		k0 = i + N*(j - 1);       % (i,j)
		k1 = i - 1 + N*(j - 1);   % (i-1,j)   
		k2 = i + 1 + N*(j - 1);   % (i+1,j)   
		k3 = i + N*(j - 2);       % (i,j-1)   
		k4 = i + N*j;             % (i,j+1)             
        k5 = i + N*(j + 1);       % (i,j+2)
        k6 = i + N*(j - 3);       % (i,j-2)
%         ka = i + N*[0:1:N-1];                 % (i,all j)
        
		if (ityp(k0) == 0)
            
			A(k0,k0) = -((2/dr2)+2/(dt2*r(i).^2));
		    A(k0,k1) = ((r(i)-dr/2)/(r(i)*dr2));
			A(k0,k2) = ((r(i)+dr/2)/(r(i)*dr2));
			A(k0,k3) = (1/(dt2*r(i).^2));
			A(k0,k4) = (1/(dt2*r(i).^2));
            
			b(k0) = 2*r(i)^2;   % right-hand side of Poisson equation
            
		  elseif (ityp(k0) == 1)
              
		     A(k0,k0) = 1;
             b(k0) = 0;
             
%              A(k0,k0)=1;                  % partc implenting different BC
%              A(k0,ka)=-1/N;
%              b(k0)=0;
             
		elseif (ityp(k0) == 2) 
            
		    A(k0,k0) = 1;
			b(k0) = fb(j);
		elseif (ityp(k0) == 3)
            
		    A(k0,k0) = -3/(2*dt);
            A(k0,k4) =  2/(dt);
            A(k0,k5) = -1/(2*dt);
            
			b(k0) = 0;
		elseif (ityp(k0) == 4)
			A(k0,k0) = 3/(2*dt);
            A(k0,k3) = -2/(dt);
            A(k0,k6) = 1/(2*dt);
            
			b(k0) = 0;
		end
	end
end

u = pinv(A)*b;

for j=1:N
	for i=1:N
		k = i + N*(j - 1);
		X1(i,j) = r(i);
		X2(i,j) = t(j);
        Y(i,j) = u(k);                           % numerical solution
		Z(i,j) = X1(i,j)^4*sin(2*X2(i,j))^2/4;   % exact solution
	end
end

X11=X1.*cos(X2);  % Convert Polar to Cartisan
X22=X1.*sin(X2);

err=abs(Y-Z);

figure;
set(gcf,'color','w');


subplot(1,3,1);
surf(X11,X22,Y);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Unumer(exactBC)','FontSize',14);
zlim([0,0.25]);

subplot(1,3,3);
surf(X11,X22,Z);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Uexact','FontSize',14);
zlim([0,0.25]);


% Numerical BC at r=0  P2.partc


for j=1:N
    
	for i=1:N
        
		k0 = i + N*(j - 1);       % (i,j)
		k1 = i - 1 + N*(j - 1);   % (i-1,j)   
		k2 = i + 1 + N*(j - 1);   % (i+1,j)   
		k3 = i + N*(j - 2);       % (i,j-1)   
		k4 = i + N*j;             % (i,j+1)             
        k5 = i + N*(j + 1);       % (i,j+2)
        k6 = i + N*(j - 3);       % (i,j-2)
        ka = i + 1 + N*[0:1:N-1];     % (i,all j)
        
		if (ityp(k0) == 0)
            
			A(k0,k0) = -((2/dr2)+2/(dt2*r(i).^2));
		    A(k0,k1) = ((r(i)-dr/2)/(r(i)*dr2));
			A(k0,k2) = ((r(i)+dr/2)/(r(i)*dr2));
			A(k0,k3) = (1/(dt2*r(i).^2));
			A(k0,k4) = (1/(dt2*r(i).^2));
            
			b(k0) = 2*r(i)^2;   % right-hand side of Poisson equation
            
		  elseif (ityp(k0) == 1)
              
% 		     A(k0,k0) = 1;
%              b(k0) = 0;
             
             A(k0,k0)=1;                  % partc implenting different BC
             A(k0,ka)=-1/(N-1);
             b(k0)=0;
             
		elseif (ityp(k0) == 2) 
            
		    A(k0,k0) = 1;
			b(k0) = fb(j);
		elseif (ityp(k0) == 3)
            
		    A(k0,k0) = -3/(2*dt);
            A(k0,k4) =  2/(dt);
            A(k0,k5) = -1/(2*dt);
            
			b(k0) = 0;
		elseif (ityp(k0) == 4)
			A(k0,k0) = 3/(2*dt);
            A(k0,k3) = -2/(dt);
            A(k0,k6) = 1/(2*dt);
            
			b(k0) = 0;
		end
	end
end

u = pinv(A)*b;

for j=1:N
	for i=1:N
		k = i + N*(j - 1);
        Y1(i,j) = u(k);                           % numerical solution
	end
end


err2=abs(Y-Z);


subplot(1,3,2);
surf(X11,X22,Y1);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Unumer (numerBC)','FontSize',14);
zlim([0,0.25]);

figure;
set(gcf,'color','w');

subplot(1,2,1)
surf(X11,X22,err);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('GTE','FontSize',14);
zlim([0,0.001]);
title('analytical BC GTE')

subplot(1,2,2)
surf(X11,X22,err2);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('GTE','FontSize',14);
zlim([0,0.001]);
title('Numerical BC GTE')


