
%%% Log normal RV : Coefficient Samples for Monte Carlo : NL Poisson

clear

mug_0 = -12.444;
sig_0 = 0.1;


mug_1 = -18.420;
sig_1 = 0.1;

n = 100000;

cd = zeros(n,2);

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 2;
dim     = 2;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))


xi1 = randn(1,n);
xi2 = randn(1,n);



psi1(:,1) = ones(1,n);
psi1(:,2) = xi1;
psi1(:,3) = (xi1.^2-1);


psi2(:,1) = ones(1,n);
psi2(:,2) = xi2;
psi2(:,3) = xi2.^2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%  Input Expansion %%%%%%%%%%%%%%%%

mul_0 = exp(mug_0 + sig_0^2/2);

ln_sd0 = sqrt(mul_0^2 * (exp(sig_0^2) -1) );

sprintf("log normal mean of RV 1 is %.8e",mul_0)
sprintf("log normal sd of of RV 1 %.8e",ln_sd0)

alpha(1) = mul_0 ;
alpha(2) = mul_0*sig_0;
alpha(3) = mul_0*sig_0^2/2;   
 

mul_1 = exp(mug_1 + sig_1^2 /2);
ln_sd1 = sqrt(mul_1^2 * (exp(sig_1^2) -1));

sprintf("log normal mean of RV 2 is %.4e",mul_1)
sprintf("log normal sd of of RV 2 is %.4e",ln_sd1)

beta(1) = mul_1;
beta(2) = mul_1*sig_1;   
beta(3) = mul_1*sig_1^2/2;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%



for j =1:n
      
       stiff1 = 0.0;
       stiff2 = 0.0;
       
       for k = 1:3
    
            stiff1 = stiff1 + alpha(k) * psi1(j,k);
            
            stiff2 = stiff2 + beta(k) * psi2(j,k);
    
       end
       
       cd(j,1) = stiff1;
       cd(j,2) = stiff2;


end

% cd = cd';

save('nu.txt','cd','-ascii');
dlmwrite('nu.txt',cd,'delimiter','\t');


