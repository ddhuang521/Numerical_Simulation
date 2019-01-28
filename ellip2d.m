
% Solve the 2D Poisson equation with Dirichlet by 2nd-order central
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ellip2d

clc; clear all; format long;

xl = 0;
xr = 1;

yl = 0;
yr = 1;

J = input('Enter number of subintervals ==> ');
h = (xr - xl)/J;

N = (J - 1)^2;

x = zeros(J-1,1);
y = zeros(J-1,1);
Uexac = zeros(J-1,J-1);

A = zeros(N,N);
F = zeros(N,1);

for i=1:N
    A(i,i) = 4;
    if (i+1 <= N) & (mod(i,J-1) ~= 0)
        A(i,i+1) = -1;
    end
    if (i+J-1 <= N)
        A(i,i+J-1) = -1;
    end
    if (i-1 >= 1) & (mod(i-1,J-1) ~= 0)
        A(i,i-1) = -1;
    end
    if (i-(J-1) >= 1)
        A(i,i-(J-1)) = -1;
    end
end

for j=1:J-1
    y(j) = j*h;
    for i= 1:J-1
        x(i) = i*h;
        F((J-1)*(j-1)+i) = h^2*(2*pi^2*sin(pi*x(i))*sin(pi*y(j)));
        Uexac(i,j) = sin(pi*x(i))*sin(pi*y(j));
    end
end

V = sparse(A)\F;

for j=1:J-1
    for i=1:J-1
        U(i,j) = V((J-1)*(j-1)+i);
    end
end

ERR = abs(U - Uexac);

subplot(1,2,1);
surf(x,y,U);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Unumer','FontSize',14);
axis([0 1 0 1 0 1]);
subplot(1,2,2);
surf(x,y,Uexac);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Uexact','FontSize',14);
axis([0 1 0 1 0 1]);

figure;

surf(x,y,ERR);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Error','FontSize',14);




