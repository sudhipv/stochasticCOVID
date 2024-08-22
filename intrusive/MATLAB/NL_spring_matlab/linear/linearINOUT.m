

%%% Linear : Monte Carlo  : INPUT OUTPUT Check


%%% ( K_0*u  ) = f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

F = 500; % kN


n = 80000000;

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 1;
dim     = 1;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% n_inpce = 2;

xi1 = randn(n,1);


%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%
mug= 100;
cov = 0.1;
sigg = cov*mug;

stiff(n,1) = 0;
stiff(:) = mug + sigg .* xi1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


U = zeros(n,1);


for j =1:n
   

        u = F/stiff(j);
        
        U(j) = u;
            

end
    


% Mean and Standard Deviation of solution
 U_mean = sum(U,1)/n
 
 sum = 0;
 
 for j = 1:n
   
     sum = sum + (U(j) - U_mean).^2;
     
     
 end
 
 U_var = sum/(n-1)
 
 U_std = sqrt(U_var)

%%% To find number of diverged samples
zz = find(~U);

cov = U_std/U_mean

figure(100)
ksdensity(stiff)

figure(200)
ksdensity(U)

