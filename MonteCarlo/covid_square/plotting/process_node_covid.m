%%% MCS COVID Model


k1 = load('./results/mcs_covid/MCS_nodeI10.dat');

k2 = load('./results/mcs_covid/MCS_nodeI11000.dat');

k3 = load('./results/mcs_covid/MCS_nodeI12000.dat');

k4 = load('./results/mcs_covid/MCS_nodeI13000.dat');

k5 = load('./results/mcs_covid/MCS_nodeI14000.dat');

% k4 = load('./output/solMCS_4.dat');

nS = 5000;

k1sum = sum(k1,2);
k2sum = sum(k2,2);
k3sum = sum(k3,2);
k4sum = sum(k4,2);
k5sum = sum(k5,2);

kfullsum = k1sum+k2sum+k3sum+k4sum+k5sum;

solmean = kfullsum/nS;

squaresum = zeros(length(k1sum),1);

for i = 1:1000
    
    squaresum = squaresum + (k1(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k2(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k3(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k4(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k5(:,i) - solmean(:,1)).^2;
    
end


solvar = squaresum/(nS-1);
solsd = sqrt(solvar);



%% 

t = linspace(0.1,1,10);


intr_node1 = load("./results/processed/intr_nodeI1.dat");
intr_node2 = load("./results/processed/intr_nodeI2.dat");

n = 300000;
xi = randn(n,1);

psi(:,1) = ones(1,n);
psi(:,2) = xi;
psi(:,3) = (xi.^2-1);
psi(:,4) = (xi.^3-3.*xi);
psi(:,5) = (xi.^4-6*xi.^2+3);

U_sample = zeros(10,n);

U_time = zeros(10,n);

for t = 1:10
    
    for pp = 1:4   
    U_sample(:,) = U_sample + intr_node1(:,pp).*psi(:,pp)/sqrt(norm_squared(pp));
    end

end

% Mean and Standard Deviation of solution
 U_mean = sum(U_nl,1)/n;
 
 sum = 0;
 
 for j = 1:n
   
     sum = sum + (U_nl(j) - U_mean).^2;
     
     
 end
 
 U_var = sum/(n-1);
 
 U_std = sqrt(U_var);





figure1 = figure;
axes1 = axes('Parent',figure1);
plot(t,solmean,'LineWidth',2);
xlabel({'t'});
ylabel({'Mean ; Infected density'});
% legend('MCS - 6000 samples');
set(axes1,'FontSize',16);
xlim([0.1,1.1])


figure2 = figure;
axes1 = axes('Parent',figure2);
plot(t,solsd,'LineWidth',2);
xlabel({'t'});
ylabel({'SD ; Infected density'});
xlim([0.1,1.1])
% legend('MCS - 6000 samples');
set(axes1,'FontSize',16);





