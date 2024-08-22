

%%% Cubic spring : Monte Carlo 

%%% Only for 1 forcing.. use _allforces for many forcing graph

%%% ( K_0*u + K_1 * u^3 ) = f

%%%% Code to check the uncertainty for ku^3 = f compared to ku = f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

F = 200; % kN
% k_0 = 50; % N/m
% k_1 = 2; % N/m

%%% Multiplier for non linear problem
mult = 0;

mug_0 = 3.91204;
mug_1 = 0.6935;

sig_0 = 0.3;
sig_1 = 0.0;

n = 20000000;

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 4;
dim     = 1;

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
psi1(:,3) = (xi1.^2-1);
psi1(:,4) = (xi1.^3-3.*xi1);
psi1(:,5) = (xi1.^4-6.*xi1.^2 + 3);


psi2(:,1) = ones(1,n);
psi2(:,2) = xi2;
psi2(:,3) = xi2.^2-1;
psi2(:,4) = xi2.^3-3.*xi2;
psi2(:,5) = (xi2.^4-6.*xi2.^2 + 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%


mul_0 = exp(mug_0 + sig_0^2/2);

ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );


sprintf("log normal mean is %.4f",mul_0)
sprintf("log normal sd is %.4f",ln_sd0)

alpha(1) = mul_0 ;
alpha(2) = mul_0*sig_0;
alpha(3) = mul_0*sig_0^2/2;   
alpha(4) = mul_0*sig_0^3/6;     
alpha(5) = mul_0*sig_0^4/24;
alpha(6) = mul_0*sig_0^5/120;
alpha(7) = mul_0*sig_0^6/720;
alpha(8) = mul_0*sig_0^7/5040;   

mul_1 = exp(mug_1 + sig_1^2 /2);
ln_sd1 = sqrt(mul_1^2 * (exp(sig_1^2) -1));

sprintf("log normal mean is %.4f",mul_1)
sprintf("log normal sd is %.4f",ln_sd1)

beta(1) = mul_1;
beta(2) = mul_1*sig_1;
beta(3) = mul_1*sig_1^2/2;   
beta(4) = mul_1*sig_1^3/6;     
beta(5) = mul_1*sig_1^4/24;
beta(6) = mul_1*sig_1^5/120;
beta(7) = mul_1*sig_1^6/720;
beta(8) = mul_1*sig_1^7/5040;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%



np = 10000;

U_nl = zeros(n,1);

iter = zeros(n,1);

inpdf(n,1) = 0;
% 
% err = zeros(size(F,2),np);

tol = 1e-6;

%%% Flag to check Picard iteration converged
flag = 0;


for j =1:n
   
       uprev = 2.0;
      
       stiff = [0,0];
       
       for k = 1:n_inpce
    
            stiff(1) = stiff(1) + alpha(k) * psi1(j,k);
%             stiff(2) = stiff(2) + beta(k) * psi1(j,k);
             inpdf(j,1) = stiff(1);
    
       end


       %%% Picard's Iteration
    for i = 1:np

        u = sqrt(F/( stiff(1)* uprev));
        
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






