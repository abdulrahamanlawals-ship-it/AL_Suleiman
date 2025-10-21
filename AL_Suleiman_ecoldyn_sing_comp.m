
% This code plots the stochastic population dynamics and singular strategy
% for sigma_e = 0.5 and omega = 13pi/20.

syms c r r0 r1 T T0 k k0 u u1 s s1 r t n n_t v r_ m w tau

k0 = 2; T = 0.5; m = 0.25; r0 = 0.1; c = 0.1; r1 = 0.1; envar = 0.5; w = 13*pi/20;

r = @(tau)r0*exp(-c/tau)+r1;
k = @(tau)k0*exp(-((tau-0.5).^2)./0.05); 

Y = @(u,s,u1,s1)(gamma(u1/s1)*gamma(1/s)*gamma((1-u1)/s1))/(gamma(u/s)*gamma(1/s1)*gamma((1-u)/s));
A = @(u,s,u1,s1)(tau^(u/s-u1/s1))*((1-tau)^((1-u)/s-((1-u1)/s1))*Y(u,s,u1,s1));

lambda = @(u,s,u1,s1)(r(tau)*(1-(A(u,s,u1,s1)*n_t/k(tau)))); % Invasion fitness
  
   tmax = 20000; tstep = 1;
  
  tau = zeros(size(tmax));
   ns = zeros(length(tmax)+1,1);
ns(1) = 2;
 
if ((T-m)>0 && (T+m)<1)  
    syms n t 

for tstep = 1:tmax 
        t = tstep;
        f = T+m*sin(w*t);
        X = envar*randn(tmax,1);                   % Coloured noise out
     sigm = 1./(1+exp(-X'));
 tau(tstep) = f*sigm(tstep);
ns(tstep+1) = ns(tstep)*exp(r(tau(tstep))*(1-(ns(tstep)/k(tau(tstep)))));

end

figure(1)
plot(ns,'k')
xlabel('Time(t)')
ylabel('Population Dynamics (N(t))')

    Tr = 101;
   n_t = ns(Tr:end-1)';
     L = tmax-100;     
    
  Time = linspace(1,tmax,tmax);
  t_ = Time(Tr:tmax);
  tau_ = tau(t_);
 
uspace = linspace(0.001,0.999,200); sspace = linspace(0.0001,2,200);
     p = length(uspace); q = length(sspace);

g1_ = zeros(p,q);                           
g2_ = zeros(p,q);

[U,S] = meshgrid(uspace,sspace); 
    
for i = 1:p
    u = uspace(i);
 
for j = 1:q
    s = sspace(j);
                                
g1 = (exp(-1./(10.*tau_))./10 + 1/10).*((n_t.*exp(20.*(tau_ - 1/2).^2).*psi(-(u - 1)./s))./(2.*s) - (n_t.*exp(20.*(tau_ - 1/2).^2).*log(1 - tau_))./(2.*s) + (n_t.*exp(20.*(tau_ - 1/2).^2).*log(tau_))./(2.*s) - (n_t.*psi(u./s).*exp(20.*(tau_ - 1./2).^2))./(2.*s));
g2 = -(exp(-1./(10.*tau_))./10 + 1/10).*((n_t.*psi(1./s).*exp(20.*(tau_ - 1/2).^2))./(2.*s.^2) - (n_t.*exp(20.*(tau_ - 1/2).^2).*log(1 - tau_).*(u - 1))./(2.*s.^2) - (n_t.*u.*psi(u./s).*exp(20.*(tau_ - 1/2).^2))./(2.*s.^2) + (n_t.*exp(20.*(tau_ - 1/2).^2).*psi(-(u - 1)./s).*(u - 1))./(2.*s.^2) + (n_t.*u.*exp(20.*(tau_ - 1/2).^2).*log(tau_))./(2.*s.^2));

g1_(i,j) = 1/L*trapz(t_,g1);
g2_(i,j) = 1/L*trapz(t_,g2);

end
end

figure(2)
[C1,h1] = contour(U,S,g1_',[0 0],'LineWidth',3);
hold on
[C2,h2] = contour(U,S,g2_',[0 0],'LineWidth',1);
xlabel('\mu (Thermal optimum)')
ylabel('\sigma (Performance Breadth)')

p1 = length(C1); q1 = length(C2);

for i1=1:p1
for j1=1:q1
   R1 = round(C1(:,i1),2,'decimals');
   R2 = round(C2(:,j1),2,'decimals'); 
if R1==R2
      mustar = R1(1);     
   sigmastar = R1(2);  
end
   
end   
end
end

Mubar = mustar
Sbar = sigmastar