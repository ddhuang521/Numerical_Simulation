clc; clear all;

k=[1,0.1,0.01]; %different epsilon to 0

for i=1:3

ur=-1;
ul=1;

y=linspace(-10,10,10000);

a=k(i);

w=ur+(ul-ur).*(1-tanh((ul-ur).*y./(4*a)))./2; 

plot(y,w);
set(gcf,'color','w');
xlabel('y','FontSize',16);
ylabel('w','FontSize',16);
hold on;

end


