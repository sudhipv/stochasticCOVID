

%%% ODE - Determinsitic case 

%%% du/dt + k(theta) u = 0

clear

k = 1;

u_0 = 1;

T = 1;
nt = 1000;
deltaT = T/nt;

t = linspace(deltaT,T,nt)';

%%% Analytical Truth

u_exact = u_0 * exp(-k.*t);


%%% Implicit Euler

u_n = u_0;

u = zeros(nt,1);

for i = 1: nt 

%%% Implicit
  u(i) = u_n/(1 + k * deltaT);
 
 %%% Explicit
%   u(i) = u_n * (1 - k * deltaT);
 
 u_n = u(i);
 

end


figure(1)
axes1 = axes('Parent',figure);
plot(t,u_exact,'LineWidth',2);
hold on
plot(t,u,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Concentration, U'});
legend('Truth','Numerical');

%%%%% Error

e = abs(u - u_exact)./u_exact;


figure(2)
axes1 = axes('Parent',figure);
plot(t,e,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Relative Error'});





