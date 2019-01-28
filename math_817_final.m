clc; clear all; format long;

% Constants
N=50;
L=16;
q=2.0;
a=0.5;
theta=1/2;

% x-space
h=L/N;
j=[-N/2:1:(N/2-1)];
x=j.*h;

% t-space
dt=0.005;
t=[0:dt:50];
nt = length(t);
r=dt/(h^2);

% matrix
S=diag(-2*ones(1,N),0)+diag(ones(1,N-1),1)+diag(ones(1,N-1),-1);
S(N,1) = 1;
S(1,N) = 1;
I=diag(ones(1,N));
U=zeros(nt,N);

% I.C.
%u = 0.5*ones(1,N);
%uf = 0.5*ones(1,N);
u = 0.5*(1+0.1*cos(pi*x/8));
%u = a*exp(1i*2*pi.*x/L);
%uf = a*exp(1i*2*pi.*x/L);


for o=1:N

   if x(o)<=8 & x(o)>=0
       
       uf(o)=0.5*(1+0.1*(1-x(o)/8));
       
   else
       
       uf(o)=0.5*(1+0.1*(1+x(o)/8));
       
   end
    
end
U(1,:) = u;
Uf(1,:) = uf;

nprint = 25;
iprint = floor(nt/nprint);

% Split Step Finite Difference

for m = 2:nt
    
    v = exp(1i*dt*q*u.*conj(u)).*u;

    w = (I-1i*r*theta*S)\((I+1i*r*(1-theta)*S)*v');
    
    u = w';
    
    U(m,:) = u;
    
%     for m=1:nt
%     
%    uc(m)=sum(abs(U(m,:)).^2*h);
%     
%      end
    
    
    	if ((m == nt) | (mod(m-1,iprint) == 0))
           plot(x,abs(U(m,:)));
           ylim([0.4,0.6]);
           title(['time t = ',num2str((m-1)*dt)]);
           
           pause(0.5);
            
            
        end

    
    
    
end

uc(1) = 0;

for m=2:nt
    
   uc(m)=abs((sum(abs(U(m,:)).^2*h)-sum(abs(U(1,:)).^2*h))/sum(abs(U(1,:)).^2*h));
    
end



figure;
surf(x,t,real(U));
zlim([-1,1]);
% 
figure;
plot(t,uc)
% 
% % Split Step Fourier Method
% % n space

% k=[-N/2:1:(N/2-1)];
% mu=2*pi.*k./L;
% nx=length(k);
% 
% 
% for m = 2:nt
%     
%     for n=1:nx
%     
%     V1(n) = (h/L)*sum(uf.*exp(-1i.*mu(n).*x));
%     U1(n) = exp(-1i*mu(n)^2*dt)*V1(n);
%    
%     end
%     
%     Un(1,:) = V1;
%     Un(m,:) = U1;
%     
%     for j=1:nx
%     
%     Uf(m,j)=sum(U1.*exp(1i.*mu.*x(j)));
%     
%     end
%     
%     uf = Uf(m,:);
%     
% end
% 
% figure;
% surf(x,t,real(Uf));
% zlim([-1,1]);
% 
% % figure;
% % plot(x,real(Uf(130,:)),x,real(U(130,:)),'ro');
% % ylim([-1,1]);
% 
% 
% figure;
% surf(k,t,real(Un));




