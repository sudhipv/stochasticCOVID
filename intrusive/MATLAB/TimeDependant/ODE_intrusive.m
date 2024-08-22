

%%% ODE - Intrusive case 

%%% du/dt + k(theta) u = 0

clear

T = 1;
nt = 1000;
deltaT = T/nt;
t = linspace(deltaT,T,nt)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%

ord_in  = 2;
ord_out = 4;
dim   = 1;
mu_g0 = 0;
sig_0 = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))

%%%%%%%%%%%G rv%%%%%%%%%%%%%%
% alpha = zeros(1,5);
% alpha(1) = 0;
% alpha(2) = 1;
%%%%%%%%%%Log NORMAL RV%%%%%%%%%%%%%%%


mul_0 = exp(mu_g0 + sig_0^2/2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_0 = zeros(n_outpce,1); % Initial condition
u_0(1) = 1;

u_n = u_0; % Previous time step

u  = zeros(n_outpce,1); % Solution at each time step
u_pc = zeros(n_outpce,nt); % Final soluton coefficient at all time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%
cijk = dlmread('./cijk_O4_D1');

ncijk = cijk(1,1);
cijk = cijk(2:size(cijk,1),:);


%%%%%%%%Coefficient MATRIX%%%%%%%%%%%%%%%%%%

Ajk = zeros(n_outpce,n_outpce);
idx = 1;


for i=1:n_inpce
    
    for j = 1:n_outpce
        
              
        for k=1:n_outpce
    
%                 for idx = 1:ncijk
                    
                    if(i == cijk(idx,1) && j == cijk(idx,2) && k == cijk(idx,3))
                                
                        Ajk(j,k) = Ajk(j,k) + cijk(idx,4) * alpha(i) * deltaT;
%                         sprintf('indices are');
%                         i
%                         j
%                         k
                        idx = idx+1; %%% use either for loop and break or  use idx +1 and no break
%                         break;
                                
                    end
                    
%                 end

    
        end
        
              
    end
    
    
end


Ajk = eye(n_outpce) +  Ajk ;



%%% Implicit Euler


for i = 1: nt 

%%% Implicit
  u = Ajk\u_n;
 
  u_n = u;
  
  u_pc(:,i) = u;
 

end

var = 0;

for k = 2:n_outpce
    
    var = var + u_pc(k,:).^2;
    
    
end

u_var = var;
u_sd = sqrt(var);


figure(1)
axes1 = axes('Parent',figure);
% plot(t,u_exact,'LineWidth',2);
% hold on
plot(t,u_pc(1,:),'LineWidth',2);
hold on
plot(t,u_pc(2,:),'LineWidth',2);
hold on
plot(t,u_pc(3,:),'LineWidth',2);
hold on
plot(t,u_pc(4,:),'LineWidth',2);
hold on
plot(t,u_pc(5,:),'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Concentration, U'});
legend('mean','1','2','3','4');





%%%%% Error

% e = abs(u_pc(1,:)' - u_exact)./u_exact;
% 
% 
figure(2)
axes1 = axes('Parent',figure);
plot(t,u_sd,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Standard Deviation'});
