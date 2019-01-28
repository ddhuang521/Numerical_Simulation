clc;clear all;

% part a

t=cputime;

N=8192;

n=[0:1:N-1];
dx=2*pi/N;
x=-pi+dx.*n;

u=exp(-x.^2);

k=[0:1:N-1];

for j=1:N
    
    ukk(j)=sum(u.*exp(-1i.*k(j).*x));
    
end

uk=abs(ukk)./N;   %sqrt(ukk.*conj(ukk))/N;

e=cputime-t

figure;
plot(x,u);
set(gcf,'color','w');
title('Gaussian Wave Packet','FontSize',16);
xlabel('x','FontSize',16);
ylabel('u','FontSize',16);
figure;
plot(k,uk);
set(gcf,'color','w');
title('DFT vs FFT','FontSize',16);
xlabel('k(0~N-1)','FontSize',16);
ylabel('|uk|','FontSize',16);
hold on;

% fft

ufft=fft(u,N);
ukfft=abs(ufft/N);

e=cputime-t

plot(k,ukfft,'ro');
legend('DFT','FFT')



