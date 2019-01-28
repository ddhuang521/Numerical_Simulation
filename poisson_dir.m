
% Solve the 2D Poisson equation with Dirichlet by 2nd-order central
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function poisson_dir

clc; clear all; format long;

L = 1;
N = input('Enter number of grid points (odd number) ==> ');
M = N*N;

x1 = linspace(0,L,N);   % x-coordinate
x2 = x1;                % y-coordinate

dx = x1(2) - x1(1);     % grid spacing
dx2 = dx*dx;

% label the grid points: boundary or interior

for j=1:N
	for i=1:N
		k = i + N*(j - 1);                    % (i,j)    
		ityp(k) = 0;                          % interior point
		if ((i == 1) & (j >= 1) & (j < N))    % left side
		    ityp(k) = 1;
		end
		if ((i == N) & (j > 1) & (j <= N))    % right side
		    ityp(k) = 2;
		end
		if ((i > 1) & (i <= N) & (j == 1))    % bottom side
		    ityp(k) = 3;
		end
		if ((i >= 1) & (i < N) & (j == N))    % top side
		    ityp(k) = 4;
		end
	end
end

% boundary conditions

fa = zeros(N,1);
fb = 0.5*x1.^2;
fc = sin(pi*x2);
fd = exp(pi)*sin(pi*x2) + 0.5*x2.^2;

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
		if (ityp(k0) == 0)
			A(k0,k0) = -4;
			A(k0,k1) = 1;
			A(k0,k2) = 1;
			A(k0,k3) = 1;
			A(k0,k4) = 1;
			b(k0) = dx2*(x1(i)^2 + x2(j)^2);   % right-hand side of Poisson equation
		elseif (ityp(k0) == 1)
		    A(k0,k0) = 1;
			b(k0) = fc(j);
		elseif (ityp(k0) == 2)
		    A(k0,k0) = 1;
			b(k0) = fd(j);
		elseif (ityp(k0) == 3)
		    A(k0,k0) = 1;
			b(k0) = fa(i);
		elseif (ityp(k0) == 4)
			A(k0,k0) = 1;
			b(k0) = fb(i);
		end
	end
end

y = sparse(A)\b;

for j=1:N
	for i=1:N
		k = i + N*(j - 1);
		X1(i,j) = x1(i);
		X2(i,j) = x2(j);
        Y(i,j) = y(k);                                                        % numerical solution
		Z(i,j) = exp(pi*X1(i,j))*sin(pi*X2(i,j)) + 0.5*(X1(i,j)*X2(i,j))^2;   % exact solution
	end
end

ERR = abs(Y - Z);

subplot(1,2,1);
surf(X1,X2,Y);
xlabel('y_1','FontSize',14);
ylabel('y_2','FontSize',14);
zlabel('Unumer','FontSize',14);
subplot(1,2,2);
surf(X1,X2,Z);
xlabel('y_1','FontSize',14);
ylabel('y_2','FontSize',14);
zlabel('Uexact','FontSize',14);

figure;

surf(X1,X2,ERR);
xlabel('y_1','FontSize',14);
ylabel('y_2','FontSize',14);
zlabel('Error','FontSize',14);




