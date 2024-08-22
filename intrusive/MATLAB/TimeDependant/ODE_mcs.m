

%%%%% Time Dependant ODE Monte Carlo


%%% du/dt + k(theta) u = 0

clear


u_0 = 1;
T = 1;
nt = 1000;
deltaT = T/nt;
t = linspace(deltaT,T,nt)';

%%% Log normal RV
ord_in  = 2;
dim =1;
mu_g0 = 0;
sig_0 = 0.1;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))

%%% MCS Samples
 
n = 40000;

xi1 = randn(1,n);


psi = zeros(n,n_inpce);
% 1-st order Hermite PC: it is always fixed to psi(1) = 1
psi(:,1) = ones;

if(ord_in > 1 || ord_in == 1 )
% 2nd order Hermite PC: it?s always fixed to psi(2) = x if (nord > 0)
psi(:,2) = xi1;

% 3rd and more order Hermite PC?s are solved by recursive formula for i = 3 : nord+1
    for i = 3:ord_in+1
        psi(:,i) = xi1' .* psi(:,i-1) - (i-2) * psi(:,i-2);
    end
end

psi1(:,1) = ones(1,n);
psi1(:,2) = xi1;
psi1(:,3) = xi1.^2-1;
psi1(:,4) = xi1.^3-3.*xi1;
psi1(:,5) = (xi1.^4-6.*xi1.^2 + 3);


%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%


mul_0 = exp(mu_g0 + sig_0^2/2);
ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );

sprintf("log normal mean is %.4f",mul_0)
sprintf("log normal sd is %.4f",ln_sd0)

alpha(1) = mul_0;
alpha(2) = mul_0*sig_0;
alpha(3) = mul_0*sig_0^2/2;   
alpha(4) = mul_0*sig_0^3/6;     
alpha(5) = mul_0*sig_0^4/24;
alpha(6) = mul_0*sig_0^5/120;
alpha(7) = mul_0*sig_0^6/720;
alpha(8) = mul_0*sig_0^7/5040;


%%% Implicit Euler

u_n = u_0;

u = zeros(nt,1);

u_MCS = zeros(nt,n);


 for j = 1:n
     
     stiff = 0;
     
     u_n = u_0;

     u = zeros(nt,1);
    
    for k = 1:n_inpce
    
            stiff = stiff + alpha(k) * psi(j,k);
            
    end

    for i = 1: nt 

    %%% Implicit
      u(i) = u_n/(1 + stiff * deltaT);

      u_n = u(i);


    end
    
    u_MCS(:,j) = u;

 end


 % Mean and Standard Deviation of solution
 U_mean = sum(u_MCS,2)/n;
 
 sum = 0;
 
 for j = 1:n
   
     sum = sum + (u_MCS(:,j) - U_mean).^2;
     
     
 end
 
 U_var = sum/(n-1);
 
 U_std = sqrt(U_var);

 
 U_CoV = U_std./U_mean;


figure(1)
axes1 = axes('Parent',figure);
plot(t,U_mean,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Concentration, U'});
legend('Truth','Numerical');
% 
% %%%%% Error
% 
% e = abs(u - u_exact)./u_exact;
% 
% 
figure(2)
axes1 = axes('Parent',figure);
plot(t,U_std,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Standard Deviation'});

U_plus2sigma = U_mean + 2 .* U_std;
U_minus2sigma = U_mean - 2 .* U_std;


figure(3)
axes1 = axes('Parent',figure);
% plot(t,U_mean,'LineWidth',2);
% hold on
% plot(t,U_plus2sigma,'LineWidth',2);
% hold on
% plot(t,U_minus2sigma,'LineWidth',2);
hold on
tn = [t' fliplr(t')];
uplusn = [U_plus2sigma', fliplr(U_minus2sigma')];
fill(tn, uplusn, 'b');
plot(t,U_mean,'LineWidth',2);

set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Standard Deviation'});


figure(4)
axes1 = axes('Parent',figure);
plot(t,U_CoV,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'CoV'});





