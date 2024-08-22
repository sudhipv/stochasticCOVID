
%%%% SIR Model - MCS

%%% Random Parameters : Log normal RV
ord_in  = 2;
dim = 1;
mu_g0 = -1.45;
sig_0 = 0.1;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))

%%% MCS Samples
 
n = 15000;

xi1 = randn(1,n);


psi1(:,1) = ones(1,n);
psi1(:,2) = xi1;
psi1(:,3) = xi1.^2-1;
psi1(:,4) = xi1.^3-3.*xi1;
psi1(:,5) = (xi1.^4-6.*xi1.^2 + 3);


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

%%% Model Parameters
gamma = 1/7;

I_0 = 0.03;
S_0 = 1- I_0;
R_0 = 0;

T = 100;
deltaT = 1/(24*2);
nt = T/deltaT;

t = linspace(deltaT,T,nt)';

tol = 1e-6;


flag = 0;

U_MCS_s = zeros(nt,n);
U_MCS_i = zeros(nt,n);
U_MCS_r = zeros(nt,n);

for j = 1:n
    
    
    U = zeros(nt,2);

    R = 0;

    U_k = [S_0;I_0];

    U_n = [S_0;I_0];

    U_0 = [S_0;I_0];

    R_n = R_0;

    U_new = zeros(2,1);
    
    beta = 0;
    
    for k = 1:n_inpce
    
            beta = beta + alpha(k) * psi1(j,k);
            
    end
    
    
    for i = 1: nt 

   
            %%% Picard Iteration
            for k = 1:50


                    U_new(1) = U_n(1)/(1 + beta * deltaT * U_k(2));

                    U_new(2) = U_n(2)/(1 + gamma* deltaT - deltaT * beta * U_new(1));

                    diff = U_new - U_k;

                    err(k,i) = (norm(diff))/norm(U_0);

                    U_k = U_new;

                    if(err(k,i) < tol)
                        iter(i) = k;
                        %sprintf("Finished picard loop with %d iteration",k)
                        flag = 1;
                        break;
                    end  

            end


            if(flag == 0)
                   sprintf("Iteration didn't converge") 
                   break;         
            end


            U_n = U_k;

            U(i,:) = U_k;

            R(i) = R_n + gamma* deltaT * U_k(2) ;

            R_n = R(i);
 
    end
    
    
    U_MCS_s(:,j) =  U(:,1);
    U_MCS_i(:,j) =  U(:,2);
    U_MCS_r(:,j) = R';
    
    
    
end

%%
% Mean and Standard Deviation of solution
 U_mean_s = sum(U_MCS_s,2)/n;
 U_mean_i = sum(U_MCS_i,2)/n;
 U_mean_r = sum(U_MCS_r,2)/n;
 
 s_s = 0;
 s_i = 0;
 s_r = 0;
 
 for j = 1:n
   
     s_s = s_s + (U_MCS_s(:,j) - U_mean_s).^2;
     s_i = s_i + (U_MCS_i(:,j) - U_mean_i).^2;
     s_r = s_r + (U_MCS_r(:,j) - U_mean_r).^2;
     
     
 end
 
 U_var_s = s_s/(n-1);
 U_var_i = s_i/(n-1);
 U_var_r = s_r/(n-1);
 
 U_std_s = sqrt(U_var_s);
 U_std_i = sqrt(U_var_i);
 U_std_r = sqrt(U_var_r);

 
 U_CoV_s = U_std_s ./U_mean_s;
 U_CoV_i = U_std_i ./U_mean_i;
 U_CoV_r = U_std_r ./U_mean_r;


figure1 = figure(1);
axes1 = axes('Parent',figure1);
plot(t,U_mean_s,t,U_mean_i,t,U_mean_r,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Proportion of mean population'});
legend('S','I','R');
grid on

% figure4 = figure(4);
% axes1 = axes('Parent',figure4);
% plot(t,U_std_s,'LineWidth',2);
% set(axes1,'FontSize',16);
% xlabel({'Time (s)'});
% ylabel({'sd'});
% legend('S');
% grid on

figure5 = figure(5);
axes1 = axes('Parent',figure5);
plot(t,U_std_s,t,U_std_i,t,U_std_r,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'sd'});
legend('S','I','R');
grid on


figure7 = figure(7);
axes1 = axes('Parent',figure7);
plot(t,U_CoV_s,t,U_CoV_i,t,U_CoV_r,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'CoV'});
legend('S','I','R');


%%%% PDF construction
%%
loc = 962;

figure8 = figure(8);
axes1 = axes('Parent',figure8);
[f_s,xi_s] = ksdensity(U_MCS_s(loc,:));
plot(xi_s,f_s, 'LineWidth',2)
hold on
[f_i,xi_i] = ksdensity(U_MCS_i(loc,:));
hold on
plot(xi_i,f_i, 'LineWidth',2)
hold on
[f_r,xi_r] = ksdensity(U_MCS_r(loc,:));
plot(xi_r,f_r, 'LineWidth',2)
set(axes1,'FontSize',16);
xlabel({'x'});
ylabel({'f(x)'});
legend('S','I','R');










