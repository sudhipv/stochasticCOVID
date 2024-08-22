

%%%% NL spring - Intrusive SSFEM using Picard's Iteration


clear

f = 200; % kN
% k_0 = 50; % N/m
% k_1 = 2; % N/m

ord_in  = 2;
ord_out = 7;
dim_i   = 1;
dim_o   = 1;


mult = 0;


mug_0 = 3.91204;
mug_1 = 0.6935;

sig_0 = 0.1;
sig_1 = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim_i)/(factorial(ord_in)*factorial(dim_i))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim_o)/(factorial(ord_out)*factorial(dim_o))




%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%

mul_0 = exp(mug_0 + sig_0^2/2);

ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );


sprintf("log normal mean is %.4f",mul_0)
sprintf("log normal sd is %.4f",ln_sd0)


alpha(1) = mul_0 ;
alpha(2) = mul_0*sig_0;
alpha(3) = mul_0*sig_0^2/sqrt(2);   
alpha(4) = mul_0*sig_0^3/sqrt(6);     
alpha(5) = mul_0*sig_0^4/sqrt(24);
alpha(6) = mul_0*sig_0^5/sqrt(120);
alpha(7) = mul_0*sig_0^6/sqrt(720);
alpha(8) = mul_0*sig_0^7/sqrt(5040);   

mul_1 = exp(mug_1 + sig_1^2 /2);
ln_sd1 = sqrt(mul_1^2 * (exp(sig_1^2) -1));

sprintf("log normal mean is %.4f",mul_1)
sprintf("log normal sd is %.4f",ln_sd1)

beta(1) = mul_1;
beta(2) = mul_1*sig_1;
beta(3) = mul_1*sig_1^2/sqrt(2);   
beta(4) = mul_1*sig_1^3/sqrt(6);     
beta(5) = mul_1*sig_1^4/sqrt(24);
beta(6) = mul_1*sig_1^5/sqrt(120);
beta(7) = mul_1*sig_1^6/sqrt(720);
beta(8) = mul_1*sig_1^7/sqrt(5040);   


%%% Gaussian RV
% 
% alpha(1) = k_0;
% alpha(2) = sig_0;
% 
% beta(1) = k_1;
% beta(2) = sig_1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%

cijk = dlmread('./../Sijklm/cijk_O10_D1');

ncijk = cijk(1,1);
cijk = cijk(2:size(cijk,1),:);

Tjklm = dlmread('./../Sijklm/Tijkl_O10_D1');
nTlm = Tjklm(1,1);
Tjklm = Tjklm(2:size(Tjklm,1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%


% uprev = ones(1,n_outpce);

% uprev = uprev .*[4,0.1,0.001,0.0001];

uprev = zeros(1,n_outpce);

% %  Assembly block struture for cijk


Clm = zeros(n_outpce);

for i = 1:n_inpce
    
    for j = 1:n_outpce
        
        for k = 1: n_outpce
        
                for idx = 1:ncijk
                    
                    if(i == cijk(idx,1) && j == cijk(idx,2) && k == cijk(idx,3))


                        Clm(j,k) = Clm(j,k) + cijk(idx,4) * alpha(i);
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

np = 50;

tol = 1e-6;

iter = 0;

diff = zeros(np,1);
% iter = zeros(size(F,2),1);
% 
% err = zeros(size(F,2),np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1 : np

     Tlm = zeros(n_outpce);

        for m = 1: n_outpce

            for l = 1:n_outpce

                for k = 1:n_outpce

                      for j = 1:n_inpce


                          for idx = 1:nTlm

                                if(j == Tjklm(idx,1) && k == Tjklm(idx,2) && l == Tjklm(idx,3) ...
                                        && m == Tjklm(idx,4))


                                    Tlm(l,m) = Tlm(l,m) + Tjklm(idx,5) * mult * beta(j) * uprev(k);
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



        stiff = Clm + Tlm;

        Unl = stiff\F';
        
        
        diff(p) = norm(uprev - Unl)/norm(Unl);
        
        
        uprev = Unl;
        
        iter = iter + 1;
       
        if(diff(p) < tol)
            break;
        end  



end

Unl

var = 0;

for k = 2:n_outpce
    
    var = var + Unl(k)^2;
    
    
end

u_var = var
u_sd = sqrt(var)
