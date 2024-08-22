

%%% NL Cubic spring : Deterministic

%%% K_0 * u + K_1 * u^3  = f

clear

F = -300 : 25 : 300; % kN %% Picard iteration doesn't converge after -475 ??

U = zeros(size(F,2),1);

k_0 = 50; % N/m

k_1 = 2; % N/m


%%% Linear Case


%%%% F Vs U - Linear Graph
U = F./k_0;

np = 1000;
%%% Non linear Case

U_nl = zeros(size(F,2),1);

iter = zeros(size(F,2),1);

err = zeros(size(F,2),np);

tol = 1e-6;


%%% Picard's Iteration
for j =1:size(F,2)

%       uprev = U(j);
      uprev = 0;
%       uprev = -30;
%       if(j < 20)
%           uprev = -10;
%       else
%           uprev = 10;
%       end

    for i = 1:np

        u = F(j)/( k_0 + k_1*uprev^2);
        
        diff = abs(abs(uprev) - abs(u))/abs(u);
        
        uprev = u;
        
        err(j,i) = diff;
        
        if(diff < tol)
            U_nl(j) = u;
            iter(j) = i;
            break;
        end  


    end


end

figure(1)
plot(U,F,'MarkerSize',15,'LineWidth',2)
hold on
plot(U_nl,F,'MarkerSize',15,'Marker','pentagram','LineWidth',2);
xlabel({'Displacement(m)'});
ylabel({'Force (kN)'});

figure(2)

E = err(j,:);

it = iter(j);

i = 1:1:it;

semilogy(i,E(1:it),'LineWidth',2);
xlabel({'Number of iterations'});
ylabel({'Error'});










