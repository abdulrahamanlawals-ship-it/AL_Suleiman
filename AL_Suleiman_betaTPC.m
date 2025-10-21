
% Illustration of thermal perf. curves for different values of mu and 
% s
 
 u1 = 0.5;  s1 = 1;  a1 = u1/s1;   b1 = (1-u1)/s1;    
 u2 = 0.95;  s2 = 0.1; a2 = u2/s2;   b2 = (1-u2)/s2;    
 u3 = 0.05;  s3 = 0.1; a3 = u3/s3;   b3 = (1-u3)/s3;    

 u4 = 0.5;  s4 = 0.1;   a4 = u4/s4;   b4 = (1-u4)/s4;    
 u5 = 0.25;  s5 = 0.075; a5 = u5/s5;   b5 = (1-u5)/s5; 
 u6 = 0.75;  s6 = 0.075; a6 = u6/s6;   b6 = (1-u6)/s6; 

T = linspace(0.001,0.999);  % Ambient temp.

% Thermal perf. curves
 f1 = betapdf(T,a1,b1);
 f2 = betapdf(T,a2,b2);   
 f3 = betapdf(T,a3,b3);  
 f4 = betapdf(T,a4,b4);
 f5 = betapdf(T,a5,b5);
 f6 = betapdf(T,a6,b6);
 
figure(1)
plot(T,f1,'k',T,f2,'r',T,f3,'b',T,f4,'g',T,f5,'b--',T,f6,'r--','LineWidth',1)
xlabel(['Ambient temperature, ',char([0xD835 0xDF0F]),''])
ylabel(['Performance, ',' \beta(' char([0xD835 0xDF0F]), ';' ,'\mu' ',' '\sigma)'])
axis([0 1 0 4])
legend('\mu = 0.5,   \it{s} = \rm{1}','\mu = 0.95, \it{s} = \rm{0.1}','\mu = 0.05, \it{s} = \rm{0.1}','\mu = 0.5,   \it{s} = \rm{0.1}','\mu = 0.25, \it{s} = \rm{0.075}','\mu = 0.75, \it{s} = \rm{0.075}','Location','northoutside')
xticks([0 0.5 1])
yticks([0 2 4])