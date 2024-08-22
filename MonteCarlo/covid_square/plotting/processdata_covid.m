%%% MCS COVID Model

clear

k1 = load('./results/mcs_covid/MCS_avgI0.dat');

k2 = load('./results/mcs_covid/MCS_avgI1000.dat');

k3 = load('./results/mcs_covid/MCS_avgI2000.dat');

k4 = load('./results/mcs_covid/MCS_avgI3000.dat');

k5 = load('./results/mcs_covid/MCS_avgI4000.dat');


nS = 5000;

nt = 10;

ksum = zeros(nt,5);
kfullsum = 0;

ksum(:,1) = sum(k1,2);
ksum(:,2) = sum(k2,2);
ksum(:,3) = sum(k3,2);
ksum(:,4) = sum(k4,2);
ksum(:,5) = sum(k5,2);

for i=1:5
    
kfullsum = kfullsum + ksum(:,i);

end

solmean = kfullsum/nS;

squaresum = zeros(nt,1);

for i = 1:1000
    
    squaresum = squaresum + (k1(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k2(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k3(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k4(:,i) - solmean(:,1)).^2;
    squaresum = squaresum + (k5(:,i) - solmean(:,1)).^2;
    
end


solvar = squaresum/(nS-1);
solsd = sqrt(solvar);



%% pdf at a point

% loc = 1300; % center

pcenter = zeros(1,nS);

pcenter(1,1:1000) = k1(loc,:);

pcenter(1,1001:2000) = k2(loc,:);

pcenter(1,2001:3000) = k3(loc,:);

% pcenter(1,4501:6000) = k4(loc,:);

%%
loc = 1800; % 2500 - 0.2,0.8, 1800 : 0.3,0.7

p2 = zeros(1,nS);

p2(1,1:1500) = k1(loc,:);

p2(1,1501:3000) = k2(loc,:);

p2(1,3001:4500) = k3(loc,:);

p2(1,4501:6000) = k4(loc,:);

% figure2 = figure;
% axes1 = axes('Parent',figure2);
ksdensity(p2);
xlabel({'x'});
ylabel({'f(x)'});
legend('MCS - 6000 samples');
set(axes1,'FontSize',16);





% k1(:,502:1000) = k2(:,1:499);

% save('./results/mcs_covid/MCSnodespaceS0.dat','k1','-ascii')
