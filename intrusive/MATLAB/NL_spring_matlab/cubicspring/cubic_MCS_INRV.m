

%%% Cubic spring : Monte Carlo 

%%% Only for 1 forcing.. use _allforces for many forcing graph

%%% ( K_0*u + K_1 * u^3 ) = f

% % % K_0 and K_1 are independant RV's

%%% Different way of input expansion, gives similar result as simple case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

F = 200; % kN
% k_0 = 50; % N/m
% k_1 = 2; % N/m

%%% Multiplier for non linear problem
mult = 1;

mug_0 = 3.91204;
mug_1 = 0.6935;

sig_0 = 0.1;
sig_1 = 0.1;

n = 20000000;

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 2;
dim     = 2;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% n_inpce = 2;

xi1 = randn(1,n);
xi2 = randn(1,n);


%%% test that both are independant normal RV
% val = xi1* xi2';
% val/n
% mean(xi1)
% mean(xi2)
% var(xi1)
% var(xi2)



psi1(:,1) = ones(1,n);
psi1(:,2) = xi1;
psi1(:,3) = xi2;
psi1(:,4) = (xi1.^2-1);
psi1(:,5) = xi1.*xi2;
psi1(:,6) = xi2.^2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%


mul_0 = exp(mug_0 + sig_0^2/2);

ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );


sprintf("log normal mean is %.4f",mul_0)
sprintf("log normal sd is %.4f",ln_sd0)

alpha(1) = mul_0 ;
alpha(2) = mul_0*sig_0;
alpha(3) = 0;   
alpha(4) = mul_0*sig_0^2/2;     
alpha(5) = 0;
alpha(6) = 0;
 

mul_1 = exp(mug_1 + sig_1^2 /2);
ln_sd1 = sqrt(mul_1^2 * (exp(sig_1^2) -1));

sprintf("log normal mean is %.4f",mul_1)
sprintf("log normal sd is %.4f",ln_sd1)

beta(1) = mul_1;
beta(2) = 0;
beta(3) = mul_1*sig_1;   
beta(4) = 0;     
beta(5) = 0;
beta(6) = mul_1*sig_1^2/2;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%



np = 10000;

U_nl = zeros(n,1);

iter = zeros(n,1);
% 
% err = zeros(size(F,2),np);

tol = 1e-6;

%%% Flag to check Picard iteration converged
flag = 0;


for j =1:n
   
       uprev = 1.0;
      
       stiff = [0,0];
       
       for k = 1:n_inpce
    
            stiff(1) = stiff(1) + alpha(k) * psi1(j,k);
            stiff(2) = stiff(2) + beta(k) * psi1(j,k);
    
       end


       %%% Picard's Iteration
    for i = 1:np

        u = F/( stiff(1) + mult * stiff(2)* uprev^2);
        
        diff = abs( abs(uprev) - abs(u) )/abs(u);
        
        uprev = u;
        
       
        if(diff < tol)
            U_nl(j) = u;
            flag = 1;
            break;
        end  


    end
    
    iter(j) = i;
    
    if(flag == 0)
        sprintf("Iteration didn't converge")
    else
        flag = 0;
    end


end



% Mean and Standard Deviation of solution
 U_mean = sum(U_nl,1)/n
 
 sum = 0;
 
 for j = 1:n
   
     sum = sum + (U_nl(j) - U_mean).^2;
     
     
 end
 
 U_var = sum/(n-1)
 
 U_std = sqrt(U_var)

%%% To find number of diverged samples
zz = find(~U_nl);






