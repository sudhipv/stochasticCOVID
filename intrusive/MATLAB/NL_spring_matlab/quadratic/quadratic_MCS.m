

%%% NL spring : Monte Carlo 

%%% ( K_0 + K_1 * u ) * u = f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

F = 200; % kN
% k_0 = 50; % N/m
% k_1 = 2; % N/m


mug_0 = 3.91204;
mug_1 = 0.6935;

sig_0 = 0.1;
sig_1 = 0.0;

mult = 0;

n = 8000000;

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 2;
dim_i     = 1;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim_i)/(factorial(ord_in)*factorial(dim_i))
% n_inpce = 2;

xi1 = randn(1,n);
xi2 = randn(1,n);


psi1(:,1) = ones(1,n);
psi1(:,2) = xi1;
psi1(:,3) = xi1.^2-1;
psi1(:,4) = xi1.^3-3.*xi1;


psi2(:,1) = ones(1,n);
psi2(:,2) = xi2;
psi2(:,3) = xi2.^2-1;
psi2(:,4) = xi2.^3-3.*xi2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%

%%% Log normal RV
mul_0 = exp(mug_0 + sig_0^2/2);

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


%%% Gaussian RV

% alpha(1) = 50;
% alpha(2) = sig_0;
% 
% beta(1) = 2;
% beta(2) = sig_1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%



np = 50;

U_nl = zeros(n,1);

% iter = zeros(size(F,2),1);
% 
% err = zeros(size(F,2),np);

tol = 1e-6;

flag = 0;


for j =1:n
   
       uprev = 0;
      
       stiff = [0,0];
       
       for k = 1:n_inpce
    
            stiff(1) = stiff(1) + alpha(k) * psi1(j,k);
            stiff(2) = stiff(2) + mult * beta(k) * psi1(j,k);
    
       end


       %%% Picard's Iteration
    for i = 1:np

        u = F/( stiff(1) + stiff(2)* uprev);
        
        diff = abs(uprev - u)/abs(u);
        
        uprev = u;
        
        
        if(diff < tol)
            U_nl(j) = u;
            flag = 1;
            break;
        end  


    end
    
    
    if(flag == 0)
        sprintf("Iteration didn't converge")
%         sprintf("Force is%d",F(f))
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









