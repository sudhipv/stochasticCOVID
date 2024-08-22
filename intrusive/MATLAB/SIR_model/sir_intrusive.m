
%%%% SIR Intrusive

clear
% beta = 0.25;
gamma = 1/7;

I_0 = 0.03;
S_0 = 1- I_0;
R_0 = 0;

T = 100;
deltaT = 1/(24*2);
nt = T/deltaT;

t = linspace(deltaT,T,nt)';

ord_in  = 2;
ord_out = 4;
dim   = 1;
mug = -1.45;
sig = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))


mul = exp(mug + sig^2 /2);
ln_sd = sqrt(mul^2 * (exp(sig^2) -1));

sprintf("log normal mean is %.4f",mul)
sprintf("log normal sd is %.4f",ln_sd)

beta(1) = mul;
beta(2) = mul*sig;
beta(3) = mul*sig^2/sqrt(2);   
beta(4) = mul*sig^3/sqrt(6);     
beta(5) = mul*sig^4/sqrt(24);
beta(6) = mul*sig^5/sqrt(120);
beta(7) = mul*sig^6/sqrt(720);
beta(8) = mul*sig^7/sqrt(5040);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%


u_0 = zeros(n_outpce,2); % Initial condition
u_0(1,1) = S_0;
u_0(1,2) = I_0;

u_n = u_0; % Previous time step

u_k = u_n; % Picard last Iteration sol

R_n = zeros(n_outpce,1);

u  = zeros(n_outpce,2); % Solution at each time step
% Final soluton coefficient at all time steps
u_s = zeros(n_outpce,nt);  % S
u_i = zeros(n_outpce,nt);  % I
u_r = zeros(n_outpce,nt);  % R

Tjklm = dlmread('./Tijkl_O10_D1');
nTlm = Tjklm(1,1);
Tjklm = Tjklm(2:size(Tjklm,1),:);

tol = 1e-6;
flag = 0;

for i = 1: nt 

    
    
    
    for p = 1 : 50
        
        
        
        Tlm_s = zeros(n_outpce);
        Tlm_i = zeros(n_outpce);

        for m = 1: n_outpce

            for l = 1:n_outpce

                for k = 1:n_outpce

                      for j = 1:n_inpce


                          for idx = 1:nTlm

                                if(j == Tjklm(idx,1) && k == Tjklm(idx,2) && l == Tjklm(idx,3) ...
                                        && m == Tjklm(idx,4))


                                    Tlm_s(l,m) = Tlm_s(l,m) + Tjklm(idx,5) * beta(j) * u_k(k,2);
                                    Tlm_i(l,m) = Tlm_i(l,m) + Tjklm(idx,5) * beta(j) * u_k(k,1);
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

        end %%%% Tlm loop finished
        
        
       Ajk_s =  eye(n_outpce) + deltaT * Tlm_s; 
       
       Ajk_i =  eye(n_outpce) + eye(n_outpce)* gamma * deltaT - deltaT * Tlm_i; 
        
        
       u(:,1) = Ajk_s\u_n(:,1);
       u(:,2) = Ajk_i\u_n(:,2);
        
       
       diff = sqrt( norm(u(:,1) - u_k(:,1))^2 + norm(u(:,2) - u_k(:,2))^2 ) ;
       den = sqrt( norm(u_0(:,1))^2 + norm(u_0(:,2))^2 );
       
       err(p,i) = diff/den;
        
        
       u_k = u;
       
       if(err(p,i) < tol)
           iter(i) = p;
%            sprintf("Finished picard loop with %d iteration",p)
           flag = 1;
           break;
       end  
       
       
        
    end %%%% Picard loop finished
    
    if(flag == 0)
           sprintf("Iteration didn't converge") 
           break;
    end
      
  %%% Implicit

 
  u_n = u;
  
  u_s(:,i) = u(:,1);
  u_i(:,i) = u(:,2);
  
  u_r(:,i) = R_n + gamma * deltaT * u(:,2);
  R_n = u_r(:,i);

end

%%

var_s = 0;
var_i = 0;
var_r = 0;

for k = 2:n_outpce
    
    var_s = var_s + u_s(k,:).^2;
    var_i = var_i + u_i(k,:).^2;
    var_r = var_r + u_r(k,:).^2;
    
    
end

u_var_s = var_s;
u_sd_s = sqrt(var_s);


u_var_i = var_i;
u_sd_i = sqrt(var_i);

u_var_r = var_r;
u_sd_r = sqrt(var_r);


figure1 = figure(1);
axes1 = axes('Parent',figure1);
plot(t,u_s(1,:),t,u_i(1,:),t,u_r(1,:),'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Proportion of mean population'});
legend('S','I','R');
grid on

% figure2 = figure(2);
% axes1 = axes('Parent',figure2);
% plot(t,u_s(2,:),t,u_i(2,:),t,u_r(2,:),'LineWidth',2);
% set(axes1,'FontSize',16);
% xlabel({'Time (s)'});
% ylabel({'PC -2'});
% legend('S','I','R');
% grid on
% 
% 
% figure3 = figure(3);
% axes1 = axes('Parent',figure3);
% plot(t,u_s(3,:),t,u_i(3,:),t,u_r(3,:),'LineWidth',2);
% set(axes1,'FontSize',16);
% xlabel({'Time (s)'});
% ylabel({'PC -3'});
% legend('S','I','R');
% grid on

% figure4 = figure(4);
% axes1 = axes('Parent',figure4);
% plot(t,u_sd_s(1,:),'LineWidth',2);
% set(axes1,'FontSize',16);
% xlabel({'Time (s)'});
% ylabel({'sd'});
% legend('S');
% grid on

figure5 = figure(5);
axes1 = axes('Parent',figure5);
plot(t,u_sd_s(1,:),t,u_sd_i(1,:),t,u_sd_r(1,:),'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'sd'});
legend('S','I','R');
grid on


%%%%%%% Mean +- 2sd S compartment
U_plus2sigma = u_s(1,:) + 2 .* u_sd_s(1,:);
U_minus2sigma = u_s(1,:) - 2 .* u_sd_s(1,:);


figure6 = figure(6);
axes1 = axes('Parent',figure6);
hold on
tn = [t' fliplr(t')];
uplusn = [U_plus2sigma, fliplr(U_minus2sigma)];
fill(tn, uplusn,[0.8 0.8 0.8]);
plot(t,u_s(1,:),'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Standard Deviation'});
legend('Mean +- 2 sigma - Susceptible')




%%%%%%% Mean +- 2sd I compartment
% U_plus2sigma = u_i(1,:) + 2 .* u_sd_i(1,:);
% U_minus2sigma = u_i(1,:) - 2 .* u_sd_i(1,:);


% figure6 = figure(6);
% axes1 = axes('Parent',figure6);
% hold on
% tn = [t' fliplr(t')];
% uplusn = [U_plus2sigma, fliplr(U_minus2sigma)];
% fill(tn, uplusn, [0.8 0.8 0.8]);
% plot(t,u_i(1,:),'LineWidth',2);
% set(axes1,'FontSize',16);
% xlabel({'Time (s)'});
% ylabel({'Standard Deviation'});
% legend('\mu +- 2\sigma - Infected')



U_CoV_s = abs( u_sd_s(1,:)./u_s(1,:) ) * 100;
U_CoV_i = abs( u_sd_i(1,:)./u_i(1,:) ) * 100;
U_CoV_r = abs ( u_sd_r(1,:)./u_r(1,:) ) * 100;


figure7 = figure(7);
axes1 = axes('Parent',figure7);
plot(t,U_CoV_s,t,U_CoV_i,t,U_CoV_r,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'CoV'});
legend('S','I','R');


%% pdf construction


lo = 3840;
n_pdf = 300000;
xi_pdf = randn(n_pdf,1);

load norm_squared070001.mat;

psi_pdf(:,1) = ones(1,n_pdf);
psi_pdf(:,2) = xi_pdf;
psi_pdf(:,3) = (xi_pdf.^2-1);
psi_pdf(:,4) = (xi_pdf.^3-3.*xi_pdf);
psi_pdf(:,5) = (xi_pdf.^4-6*xi_pdf.^2+3);


U_pdf_s = zeros(n_pdf,1);  
U_pdf_i = zeros(n_pdf,1); 
U_pdf_r = zeros(n_pdf,1);  

for pp = 1:4   
U_pdf_s(:,1) = U_pdf_s + u_s(pp,lo).*psi_pdf(:,pp)/sqrt(norm_squared(pp));
U_pdf_i(:,1) = U_pdf_i + u_i(pp,lo).*psi_pdf(:,pp)/sqrt(norm_squared(pp));
U_pdf_r(:,1) = U_pdf_r + u_r(pp,lo).*psi_pdf(:,pp)/sqrt(norm_squared(pp));
end

[f_s,xi_s] = ksdensity(U_pdf_s);
[f_i,xi_i] = ksdensity(U_pdf_i);
[f_r,xi_r] = ksdensity(U_pdf_r);

figure1 = figure();
axes1 = axes('Parent',figure1);
% plot(xi_s,f_s, 'LineWidth',2)
% hold on
% plot(xi_i,f_i, 'LineWidth',2)
% hold on 
plot(xi_r,f_r, 'LineWidth',2)
xlabel({'x'});
ylabel({'f(x)'});
legend('Intrusive - 4^{th} order');
set(axes1,'FontSize',16);








