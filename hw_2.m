
% Solve the 2D Poisson equation with Dirichlet by 2nd-order central
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; format long;

L = 1;
N = 51;
M = N*N;
J = N+1;

x = linspace(0,L,N);   % x-coordinate
y = linspace(0,pi/2,N);                % y-coordinate

dx = x(2) - x(1);     % grid spacing
dx2 = dx*dx;
dy = y(2) - y(1);
dy2 = dy*dy;

% label the grid points: boundary or interior

for j=1:N
	for i=1:N
		k = i + N*(j - 1);                    % (i,j)    
		ityp(k) = 0;                          % interior point
		if ((i == 1) & (j >= 1) & (j < N))    % left side
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
fb = sin(2.*y).^2./4;
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
		if (ityp(k0) == 0)
			A(k0,k0) = -((2/dx2)+2/(dy2*x(i).^2));
		    A(k0,k1) = ((x(i)-dx/2)/(x(i)*dx2));
			A(k0,k2) = ((x(i)+dx/2)/(x(i)*dx2));
			A(k0,k3) = (1/(dy2*x(i).^2));
			A(k0,k4) = (1/(dy2*x(i).^2));
            
			b(k0) = 2*x(i)^2;   % right-hand side of Poisson equation
		% elseif (ityp(k0) == 1)
		    %A(k0,k0) = 1;
			%b(k0) = fc(j);
		elseif (ityp(k0) == 2) 
		    A(k0,k0) = 1;
			b(k0) = fb(j);
		elseif (ityp(k0) == 3)
		    A(k0,k0) = -3/(2*dy);
            A(k0,k4) =  2/(dy);
            A(k0,k5) = -1/(2*dy);
            
			b(k0) = 0;
		elseif (ityp(k0) == 4)
			A(k0,k0) = 3/(2*dy);
            A(k0,k3) = -2/(dy);
            A(k0,k6) = 1/(2*dy);
            
			b(k0) = 0;
		end
	end
end

u = pinv(A)*b;

for j=1:N
	for i=1:N
		k = i + N*(j - 1);
		X1(i,j) = x(i);
		X2(i,j) = y(j);
        Y(i,j) = u(k);                         % numerical solution
		Z(i,j) = X1(i,j)^4*sin(2*X2(i,j))^2/4;   % exact solution
	end
end

ERR = abs(Y - Z);

%subplot(1,2,1);
%surf(X1,X2,Y);
%xlabel('y_1','FontSize',14);
%ylabel('y_2','FontSize',14);
%zlabel('Unumer','FontSize',14);
%subplot(1,2,2);
%surf(X1,X2,Z);
%xlabel('y_1','FontSize',14);
%ylabel('y_2','FontSize',14);
%zlabel('Uexact','FontSize',14);

%figure;

surf(X1,X2,ERR);
xlabel('y_1','FontSize',14);
ylabel('y_2','FontSize',14);
zlabel('Error','FontSize',14);




