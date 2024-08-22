
%%%% SIR model - deterministic - Implicit Method - Verified

beta = 0.2357;
gamma = 1/7;

I_0 = 0.03;
S_0 = 1- I_0;
R_0 = 0;

T = 100;
deltaT = 1/(24*2);
nt = T/deltaT;

t = linspace(deltaT,T,nt)';

%%% Implicit Euler


U = zeros(nt,2);

R = 0;

U_k = [S_0;I_0];

U_n = [S_0;I_0];

U_0 = [S_0;I_0];

R_n = R_0;

U_new = zeros(2,1);

tol = 1e-6;

flag = 0;

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
%             sprintf("Finished picard loop with %d iteration",k)
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

sprintf("Reproduction number is %d", beta/gamma);

figure(1)
axes1 = axes('Parent',figure);
plot(t,U(:,1),'LineWidth',2);
hold on
plot(t,U(:,2),'LineWidth',2);
hold on
plot(t,R,'LineWidth',2);
set(axes1,'FontSize',16);
xlabel({'Time (s)'});
ylabel({'Proportion of population'});
grid on

%%
%%%% Verification from colab
% % sus_colab = load('sus_colab.txt');
% % inf_colab = load('inf_colab.txt');
% % rem_colab = load('rem_colab.txt');
% % hold on
% % plot(t,sus_colab,'k--','LineWidth',2);
% % hold on
% % plot(t,inf_colab,'k--','LineWidth',2);
% % hold on
% % plot(t,rem_colab,'k--','LineWidth',2);
% % 
% % legend('Susceptible','Infected', 'Removed','Sus-Colab','Inf-Colab', 'Rem-Colab');



