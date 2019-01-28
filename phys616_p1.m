clc;clear all;

x=linspace(0,1,1000);

r=1.32e-3;
k1=[-3:1:3];
k2=[-5:1:5];
k3=[-500:1:500];

Z1=0;
Z2=0;
Z3=0;
u1=0;
u2=0;
u3=0;

for i=1:length(k1)
    
    Z1=Z1+exp(k1(i).*x);

end

for i=1:length(k2)
    
    Z2=Z2+exp(k2(i).*x);

end

for i=1:length(k3)
    
    Z3=Z3+exp(k3(i).*x);

end


for i=1:length(k1)

    u1=u1+k1(i).*exp(k1(i).*x)./Z1;

end

for i=1:length(k2)

    u2=u2+k2(i).*exp(k2(i).*x)./Z2;

end

for i=1:length(k3)

    u3=u3+k3(i).*exp(k3(i).*x)./Z3;

end

% plot(x,u1);
% hold on;
% plot(x,u2);
% hold on;
plot(x,u3);
hold on;

% exact solution

x=linspace(3,5,100);

ke=[-1:1:1];
kn=[-500:1:500];

Ze=0;
Zn=0;
ue=0;
un=0;

for i=1:length(ke)
    
    Ze=Ze+exp(ke(i).*x);

end

for i=1:length(ke)

    ue=500*(coth(x)-1./x);

end

for i=1:length(kn)
    
    Zn=Zn+exp(kn(i).*x);

end

for i=1:length(kn)

    un=un+kn(i).*exp(kn(i).*x)./Zn;

end


plot(x,ue,'o');
hold on;
%plot(x,un);
%hold on;



