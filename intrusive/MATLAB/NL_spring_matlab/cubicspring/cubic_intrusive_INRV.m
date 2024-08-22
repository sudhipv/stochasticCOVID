

%%%% NL spring - Intrusive SSFEM using Picard's Iteration

%%% Using 2 INdependant Random variables, by expanding the input with 2 RV 
%%%% exapansion and filling out coefficients corresponding to only the
%%%% respective terms.

%%% Working tested for two cases 0.1 and 0.2

clear

f = 200; % kN
% k_0 = 50; % N/m
% k_1 = 2; % N/m

% U_LDT = f/k_0;

% U_NLDT = f/(k_0 + k_1);

ord_in  = 2;
ord_out = 3;
dim   = 2;

%%% Multiplier for nonlinearity, 0 means linear
mult = 1;


mu_g0 = 3.91204;
sig_0 = 0.2;

mu_g1 = 0.6935;
sig_1 = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))




%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%

mul_0 = exp(mu_g0 + sig_0^2/2);
ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );

sprintf("log normal mean is %.4f",mul_0)
sprintf("log normal sd is %.4f",ln_sd0)

alpha(1) = mul_0 ;
alpha(2) = mul_0*sig_0;
alpha(3) = 0;   
alpha(4) = mul_0*sig_0^2/2;     
alpha(5) = 0;
alpha(6) = 0;

mul_1 = exp(mu_g1 + sig_1^2/2);
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

cijk = dlmread('./../Sijklm/cijk_O3_D2');

ncijk = cijk(1,1);
cijk = cijk(2:size(cijk,1),:);

Sjklmn = dlmread('./../Sijklm/Sijklm_O3_D2');
nSkn = Sjklmn(1,1);
Sjklmn = Sjklmn(2:size(Sjklmn,1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uprev = ones(1,n_outpce);

uprev = zeros(1,n_outpce);
uprev(1,1) = 3;

% %  Assembly block struture for cijk

Ckn = zeros(n_outpce);

for i = 1:n_inpce
    
    for j = 1:n_outpce
        
        for k = 1: n_outpce
        
                for idx = 1:ncijk
                    
                    if(i == cijk(idx,1) && j == cijk(idx,2) && k == cijk(idx,3))


                        Ckn(j,k) = Ckn(j,k) + cijk(idx,4) * alpha(i);
%                         sprintf("indices are");
%                         i
%                         j
%                         k
                        break;
                    end

                 end
            
            
        end
        
    end
    
end



Unl = zeros(1, n_outpce);

F = zeros(1,n_outpce);

F(1,1) = f;

np = 10000;

tol = 1e-6;

iter = 0;

flag = 0;

diff = zeros(np,1);
% iter = zeros(size(F,2),1);
% 
% err = zeros(size(F,2),np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1 : np


     Skn = zeros(n_outpce);

     for n = 1:n_outpce
         
        for k = 1: n_outpce

            for l = 1:n_outpce

                for m = 1:n_outpce

                      for j = 1:n_inpce


                          for idx = 1:nSkn

                                if(j == Sjklmn(idx,1) && k == Sjklmn(idx,2) && l == Sjklmn(idx,3) ...
                                        && m == Sjklmn(idx,4) && n == Sjklmn(idx,5))


                                    Skn(m,n) = Skn(m,n) + Sjklmn(idx,6) * mult * beta(j) * uprev(k) * uprev(l);
%                                     sprintf("indices are");
%                                     j
%                                     k
%                                      l
%                                      m
%                                      beta(j)
                                    break;

                                end

                          end



                      end

                end

            end

        end
        
     end



        stiff = Ckn + Skn;

        Unl = F/stiff;
        
        
        diff(p) = norm(uprev - Unl)/norm(Unl);
        
        
        uprev = Unl;
        
        iter = iter + 1;
       
        if(diff(p) < tol)
            flag = 1;
            break;    
        end  


end

if(flag == 0)
            sprintf("Iteration didn't converge")
end

Unl

var = 0;

for k = 2:n_outpce
    
    var = var + Unl(k)^2;
    
    
end

u_var = var
u_sd = sqrt(var)
