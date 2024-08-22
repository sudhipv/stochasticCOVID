

k1 = readtable('./output/u0_intr.csv');

k2 = readtable('./output/u1_intr.csv');

k3 = readtable('./output/u2_intr.csv');

k4 = readtable('./output/u3_intr.csv');


%%

u(:,1) = k1.scalaire1;
u(:,2) = k2.scalaire1;
u(:,3) = k3.scalaire1;
u(:,4) = k4.scalaire1;



lo = 1800;
n_pdf = 300000;
xi_pdf = randn(n_pdf,1);

psi_pdf(:,1) = ones(1,n_pdf);
psi_pdf(:,2) = xi_pdf;
psi_pdf(:,3) = (xi_pdf.^2-1);
psi_pdf(:,4) = (xi_pdf.^3-3.*xi_pdf);
psi_pdf(:,5) = (xi_pdf.^4-6*xi_pdf.^2+3);

load norm_squared030001.mat
U_pdf = zeros(n_pdf,1);  

for pp = 1:4   
U_pdf(:,1) = U_pdf + u(lo,pp).*psi_pdf(:,pp)/sqrt(norm_squared(pp));
end

[f_IN,xi_IN] = ksdensity(U_pdf);

figure1 = figure();
axes1 = axes('Parent',figure1);
plot(xi_IN,f_IN, 'LineWidth',2)
xlabel({'x'});
ylabel({'f(x)'});
legend('Intrusive - 3^{rd} order');
set(axes1,'FontSize',16);




