

%%% Cubic spring : Monte Carlo 

%%% ( K_0*u + K_1 * u^3 ) = f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

F = 20:20:200; % kN
% F(11) = 0.1;

% k_0 = 50; % N/m
% k_1 = 2; % N/m

% U = zeros(size(F,2),1);
% U = F./k_0;
% 
% UN = zeros(size(F,2),1);
% UN = F./(k_0+k_1);

mug_0 = 3.91204;
mug_1 = 0.6935;

sig_0 = 0.1;
sig_1 = 0.1;


%%%% Multiplier to nonlinear term, 0 means linear problem
mult = 0; 

n = 8000000;

%%%% MCS samples %%%%%%%%%%%%

ord_in  = 2;
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
psi1(:,3) = xi1.^2-1;
psi1(:,4) = xi1.^3-3.*xi1;


psi2(:,1) = ones(1,n);
psi2(:,2) = xi2;
psi2(:,3) = xi2.^2-1;
psi2(:,4) = xi2.^3-3.*xi2;

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



np = 500;

U_nl = zeros(n,size(F,2));

iter = zeros(n,1);
 
err = zeros(np,1);

tol = 1e-6;

%%% Flag to check Picard iteration converged
flag = 0;

for f=1:size(F,2)

for j =1:n
   
       uprev = 0;
      
       stiff = [0,0];
       
       for k = 1:n_inpce
    
            stiff(1) = stiff(1) + alpha(k) * psi1(j,k);
            
%%% Two independant RVs
            stiff(2) = stiff(2) + beta(k) * psi2(j,k);
            
%%% Single RV case
%% 
%             stiff(2) = stiff(2) +  beta(k) * psi1(j,k);
    
       end


       %%% Picard's Iteration
    for i = 1:np

        u = F(f)/( stiff(1) + mult * stiff(2)* uprev^2);
        
        diff = abs( abs(uprev) - abs(u) )/abs(u);
        
        uprev = u;
        
        err(i) = diff;
        
        if(diff < tol)
            U_nl(j,f) = u;
            flag = 1;
            break;
        end  


    end
    
    iter(j) = i;
    
    if(flag == 0)
        sprintf("Iteration didn't converge")
        sprintf("Force is%d",F(f))
    else
        flag = 0;
    end
    


end



% Mean and Standard Deviation of solution
 U_mean(f) = sum(U_nl(:,f))/n;
 
 sumsol = 0;
 
 for j = 1:n
   
     sumsol = sumsol + (U_nl(j,f) - U_mean(f)).^2;
     
     
 end
 
 U_var(f) = sumsol/(n-1);
 
 U_std(f) = sqrt(U_var(f));
 
 U_cov(f) = U_std(f)/U_mean(f);



end

U_mean(f)
U_std(f)
U_cov(f)



figure(1)
axes1 = axes('Parent',figure);
plot(U_mean,F,'MarkerSize',15,'Marker','pentagram','LineWidth',2);
xlabel({'Displacement(m)'});
ylabel({'Force (kN)'});
set(axes1,'FontSize',16);
legend('Linear Spring,V_0 = 0.1,1RV');

figure(2)
axes1 = axes('Parent',figure);
plot(U_std,F,'MarkerSize',15,'Marker','pentagram','LineWidth',2);
xlabel({'Standard Deviation'});
ylabel({'Force (kN)'});
set(axes1,'FontSize',16);
legend('Linear Spring,V_0 = 0.1, 1RV');

figure(3)
axes1 = axes('Parent',figure);
plot(U_cov,F,'MarkerSize',15,'Marker','pentagram','LineWidth',2);
xlabel({'CoV - Output'});
ylabel({'Force (kN)'});
set(axes1,'FontSize',16);
legend('Linear Spring,V_0 = 0.1, 1RV');


