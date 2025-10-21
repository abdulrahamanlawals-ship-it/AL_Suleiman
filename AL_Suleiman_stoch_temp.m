
% This code plots different stochastic temperature dynamics following white noise 

T0 = 0.5; m = 0.25; w = 13*pi/20; envar = 1; tmax = 100; 

t = 1:1:tmax;
f = T0+m*sin(w*t);

X = envar.*randn(length(t),1); % White noise
sigm = 1./(1+exp(-X'));

Temp = f.*sigm;

figure(1)
plot(t,Temp,'r','LineWidth',1.5)
set(gca,'fontsize',12);
xlim([0 tmax])
xticks([0 50 100])
ylim([0 0.6])
yticks([0 0.3 0.6])