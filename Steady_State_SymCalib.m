
function [SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new] = Steady_State_SymCalib(x)

global beta delta theta eta n sigma omega_c omega_k omega alpha theta_s fe fp f0 f1 xi_0 xi_1 tau S  L ...
     kappa lambda sigma_m Exporter_ROW New_exporter_ROW EX_new_ROW Exporter_Prem_ROW TO_ROW A

index = 1;

C = x(index); index = index+1;
W = x(index); index = index+1;
NT = x(index); index = index+1;

%Variables
P = zeros(1,S);
K = zeros(1,S);
P_K = zeros(1,S);
R_K = zeros(1,S);
X_C = zeros(1,S); 
I  = zeros(1,S);
D_K  = zeros(1,S);
X_K = zeros(1,S,S);
X = zeros(1,S);
D = zeros(1,S);
MC = zeros(1,S);
PI_0 = zeros(1,S);
a_p = zeros(1,S);
a_0 = zeros(1,S);
a_1 = zeros(1,S);
Psi_p = zeros(1,S);
Psi_0 = zeros(1,S);
Psi_1 = zeros(1,S);
n_p = zeros(1,S);
n_0 = zeros(1,S);
n_1 = zeros(1,S);
EV_0 = zeros(1,S); 
EV_1 = zeros(1,S); 
N = zeros(1,S); 
NE = zeros(1,S); 
NX = zeros(1,S); 
NX_1 = zeros(1,S); 
NX_0 = zeros(1,S); 
N_0 = zeros(1,S); 
Lp = zeros(1,S); 
Lf = zeros(1,S); 
TauR = zeros(1,S); 
EX = zeros(1,S); 
PI = zeros(1,S); 
EV_inf = zeros(1,S); 
P_M = zeros(1,S);
m = zeros(1,S);

X_M = zeros(1,S,S);
M = zeros(1,S);
D_M = zeros(1,S);
MC_CD = zeros(1,S);
CD = zeros(1,S);

P(1,1) = x(index); index = index+1;
K(1,1) = x(index); index = index+1;
a_0(1,1) = x(index); index = index+1;
M(1,1) = x(index); index = index+1;

for i = 2:S
    P(1,i) = P(1,1);
    K(1,i) = K(1,1);
    a_0(1,i) = a_0(1,1);
    M(1,i) = M(1,1);
end

fp(1,1) = x(index); index = index+1;
f0(1,1) = x(index); index = index+1;
f1(1,1) = f0(1,1); 
xi_0(1,1,2) = x(index); index = index+1;
xi_1(1,1,2) = xi_0(1,1,2);

for i = 2:S
fp(1,i) = fp(1,1);
f0(1,i) = f0(1,1);
f1(1,i) = f1(1,1);
xi_0(1,i,2) = xi_0(1,1,2); 
xi_1(1,i,2) = xi_1(1,1,2); 
end

fp(2,:) = fp(1,:);
f0(2,:) = f0(1,:);
f1(2,:) = f1(1,:);
xi_0(2,:,1) = xi_0(1,:,2); 
xi_1(2,:,1) = xi_1(1,:,2); 


%Consumption price
P_C = 1;
%Nominal expenditure in C
D_C = P_C*C;

%Stochastic discount factor
Q = beta;

%Price of sector-materials
for i = 1:S
for j = 1:S
    P_M(1,i) = P_M(1,i) + lambda(1,i,j)*P(1,j)^(1-sigma_m(1,i));
end
    P_M(1,i) = P_M(1,i)^(1/(1-sigma_m(1,i)));
end

%Price of sector-investment
for i = 1:S
for j = 1:S
    P_K(1,i) = P_K(1,i) + omega_k(1,i,j)*P(1,j)^(1-sigma(1,i));
end
    P_K(1,i) = P_K(1,i)^(1/(1-sigma(1,i)));
end

for i=1:S
%Interest rate
R_K(1,i) = (1/beta - (1-delta(1,i)))*P_K(1,i);
%Demand for sector-intermediates: consumption
X_C(1,i) = omega_c(1,i)*P(1,i)^(-theta)*P_C^(theta-1)*D_C;
%Investment
I(1,i) = K(1,i)-(1-delta(1,i))*K(1,i);  
%Nominal expenditure on investment
D_K(1,i) = I(1,i)*P_K(1,i);
%Nominal expenditure on materials
D_M(1,i) = P_M(1,i)*M(1,i);
end

%Demand for sector-investment
for i = 1:S
for j = 1:S
    X_K(1,i,j) = omega_k(1,i,j)*P(1,j)^(-sigma(1,i))*P_K(1,i)^(sigma(1,i)-1)*D_K(1,i);
end
end

%Demand for sector-intermediates: sector-materials
for i = 1:S
for j = 1:S
    X_M(1,i,j) = lambda(1,i,j)*P(1,j)^(-sigma_m(1,i))*P_M(1,i)^(sigma_m(1,i)-1)*D_M(1,i);
end
end

%Total demand for sector intermediates
for i = 1:S
X(1,i) = X_C(1,i);
for j = 1:S
   X(1,i) = X(1,i) + X_K(1,j,i) + X_M(1,j,i);
end
end

for i = 1:S
%Nominal expenditure on sector intermediates
D(1,i) = P(1,i)*X(1,i);

%Heterogeneous firms
%Marginal cost
MC_CD(1,i) = 1/A(1,i)*(R_K(1,i)/alpha(1,i))^alpha(1,i)*(W/kappa(1,i))^kappa(1,i)*(P_M(1,i)/(1-alpha(1,i)-kappa(1,i)))^(1-alpha(1,i)-kappa(1,i));
MC(1,i) = MC_CD(1,i);
end

for i = 1:S

%Profits constant
PI_0(1,i) = 1/theta_s(1,i)*(theta_s(1,i)*MC(1,i)/(theta_s(1,i)-1))^(1-theta_s(1,i));
%Marginal producer: productivity
a_p(1,i) = W*fp(1,i)/(PI_0(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1));
%Marginal continuing exporter
a_1(1,i) = a_0(1,i)*(xi_0(1,i,2)/xi_1(1,i,2))^(1-theta_s(1,i)) - ...
           W*(f0(1,i)-f1(1,i))*tau(1,i,2)^theta_s(1,i)*P(1,i)^(1-theta_s(1,i))/(PI_0(1,i)*D(1,i)*(1-omega(1,i,1))*(xi_1(1,i,2))^(1-theta_s(1,i)));

%Distributions
Psi_p(1,i) = eta/(eta-1)*a_p(1,i)^(1-eta);
Psi_0(1,i) = eta/(eta-1)*a_0(1,i)^(1-eta);
Psi_1(1,i) = eta/(eta-1)*a_1(1,i)^(1-eta);

n_p(1,i) = a_p(1,i)^(-eta);
n_0(1,i) = a_0(1,i)^(-eta);
n_1(1,i) = a_1(1,i)^(-eta);

%Expected values
VV = [1,  -n(1,1)*Q*n_1(1,i), -n(1,1)*Q*(1-n_1(1,i));...
      0, 1-n(1,1)*Q*n_1(1,i), -n(1,1)*Q*(1-n_1(1,i)); ...
      -n(1,1)*Q*n_0(1,i), 0, 1-(1-n_0(1,i))*n(1,1)*Q];

EE = [PI_0(1,i)*(P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1)*Psi_p(1,i)+xi_0(1,i,2)^(1-theta_s(1,i))*tau(1,i,2)^(-theta_s(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*(1-omega(1,i,1))*Psi_0(1,i))-W*fp(1,i)*n_p(1,i)-W*f0(1,i)*n_0(1,i); 
      PI_0(1,i)*(P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1)*Psi_p(1,i)+xi_1(1,i,2)^(1-theta_s(1,i))*tau(1,i,2)^(-theta_s(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*(1-omega(1,i,1))*Psi_1(1,i))-W*fp(1,i)*n_p(1,i)-W*f1(1,i)*n_1(1,i);
      PI_0(1,i)*(P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1)*Psi_p(1,i))-W*fp(1,i)*n_p(1,i)];

EVV = VV^(-1)*EE;
EV_0(1,i) = EVV(1,1);
EV_1(1,i) = EVV(2,1);
EV_inf(1,i) = EVV(3,1);

%Masses of firms

N(1,i) = NT/S;
NE(1,i) = (1-n(1,1))*N(1,i);

AA = [1, -1, -1, 0; ...
      -n(1,1)*n_1(1,i), 1, 0, 0; ...
      0, 0, 1, -n(1,1)*n_0(1,i); ...
      1, 0, 0, 1];

NN = [0; 0; n_0(1,i)*NE(1,i); N(1,i)];
NXS = AA^(-1)*NN;

%Exporters
NX(1,i)   = NXS(1,1);
%Continuation exporters
NX_1(1,i) = NXS(2,1);
%New exporters
NX_0(1,i) = NXS(3,1);
%Non exporters
N_0(1,i)  = NXS(4,1);

%Production labor
Lp(1,i) = (theta_s(1,i)-1)*PI_0(1,i)*(kappa(1,i))/W*(N(1,i)*Psi_p(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1) + ...
           tau(1,i,2)^(-theta_s(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*(1-omega(1,i,1))*(xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i)+ ...
           xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i)));

%Labor for fixed costs
Lf(1,i) = N(1,i)*n_p(1,i)*fp(1,i) + f0(1,i)*NX_0(1,i) + f1(1,i)*NX_1(1,i) + NE(1,i)*fe(1,1);

%Tariff revenue
TauR(1,i) = (1-omega(1,i,1))*(tau(1,i,2)-1)*theta_s(1,i)*PI_0(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*tau(1,i,2)^(-theta_s(1,i))*...
            (xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i)+ ...
             xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i));

%Exports
EX(1,i) =   (1-omega(1,i,1))*theta_s(1,i)*PI_0(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*tau(1,i,2)^(-theta_s(1,i))*...
            (xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i)+...
             xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i));

%Profits
PI(1,i) = PI_0(1,i)*(N(1,i)*Psi_p(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1) + ...
          tau(1,i,2)^(-theta_s(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*(1-omega(1,i,1))*...
          (xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i)+...
           xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i)))- W*Lf(1,i);
       
end

for i = 1:S
m(1,i) = (1-alpha(1,i)-kappa(1,i))/kappa(1,i)*W/P_M(1,i)*Lp(1,i);
CD(1,i) = K(1,i)^alpha(1,i)*Lp(1,i)^kappa(1,i)*m(1,i)^(1-alpha(1,i)-kappa(1,i));
PI_CD(1,i) = MC_CD(1,i)*CD(1,i) - R_K(1,i)*K(1,i) - W*Lp(1,i) - P_M(1,i)*m(1,i);

end

%%Totals
TauT = sum(TauR(1,:));
EX_T = sum(EX(1,:));
PI_T = sum(PI(1,:));
PI_CD_T = sum(PI_CD(1,:));

%% Check Walras: Consumer's Budget Constraint

WL = P_C(1,1)*C - TauT - PI_T - W*L(1,1);

for i= 1:S
   WL = WL + P_K(1,i)*I(1,i) - R_K(1,i)*K(1,i);    
end


%% Solve
index = 1;

%Free entry
SS0(index) = W*fe(1,1);
for i=1:S
    SS0(index) = SS0(index) - Q*EV_inf(1,i)/S;
end
index = index + 1;

%Labor equilibrium
SS0(index) = L(1,1) - sum(Lp(1,:)) - sum(Lf(1,:)); index = index + 1;

%Consumption prices
P_Ceq = 0;
for i = 1:S
   P_Ceq = P_Ceq + omega_c(1,i)*P(1,i)^(1-theta);
end
P_Ceq = P_Ceq^(1/(1-theta));

SS0(index) = P_C - P_Ceq; index = index+1;

for i = 1:1
    %Marginal new exporter
    SS0(index) = W*f0(1,i) - PI_0(1,i)*a_0(1,i)*xi_0(1,i,2)^(1-theta_s(1,i))*D(1,i)*(1-omega(1,i,1))/(tau(1,i,2)^theta_s(1,i)*P(1,i)^(1-theta_s(1,i))) - ...
            n(1,1)*Q*(EV_1(1,i)-EV_inf(1,i)); index = index+1;
    %Sector prices
   SS0(index) = P(1,i)^(1-theta_s(1,i)) - (omega(1,i,1)*N(1,i)*(theta_s(1,i)*MC(1,i)/(theta_s(1,i)-1))^(1-theta_s(1,i))*Psi_p(1,i) + ...
       (1-omega(1,i,1))*(tau(1,i,2)*theta_s(1,i)/(theta_s(1,i)-1)*MC(1,i))^(1-theta_s(1,i))*...
       (xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i)+ ...
        xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i))); index = index + 1;
    %Capital equilibrium
    SS0(index) = K(1,i) - alpha(1,i)/kappa(1,i)*W/R_K(1,i)*Lp(1,i); index = index+1; 
    %Manufacturing
    SS0(index) = M(1,i) - m(1,i); index = index+1;
end

%%Additional moments
for i = 1:S
     EX_1(1,i) =   (1-omega(1,i,1))*theta_s(1,i)*PI_0(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*tau(1,i,2)^(-theta_s(1,i))*...
            (xi_1(1,i,2)^(1-theta_s(1,i))*NX_1(1,i)/n_1(1,i)*Psi_1(1,i));
     EX_0(1,i) =   (1-omega(1,i,1))*theta_s(1,i)*PI_0(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*tau(1,i,2)^(-theta_s(1,i))*...
            (xi_0(1,i,2)^(1-theta_s(1,i))*NX_0(1,i)/n_0(1,i)*Psi_0(1,i));
     GO(1,i)    = theta_s(1,i)*PI_0(1,i)*N(1,i)*Psi_p(1,i)*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1) + EX(1,i);
     GO_EX(1,i) = theta_s(1,i)*PI_0(1,i)*(NX_1(1,i)*Psi_1(1,i)+NX_0(1,i)*Psi_0(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1) + EX(1,i);
     GO_NoEX(1,i) = theta_s(1,i)*PI_0(1,i)*(N(1,i)*Psi_p(1,i)-NX_1(1,i)*Psi_1(1,i)-NX_0(1,i)*Psi_0(1,i))*P(1,i)^(theta_s(1,i)-1)*D(1,i)*omega(1,i,1);
     EX_GO(1,i) = EX(1,i)/GO(1,i);
     TO(1,i) = 2*(EX(1,i))/GO(1,i);
     Exporter(1,i) = 1 - N_0(1,i)/N(1,i);
     New_exp(1,i) = NX_0(1,i)/NX(1,i);
     Exporter_Prem(1,i) = GO_EX(1,i)/GO(1,i)/Exporter(1,i); 
     EX_new(1,i) = EX_0(1,i)/EX(1,i);
end

for i = 1:1
    SS0(index) = Exporter(1,i) - Exporter_ROW(1,i); index = index+1;
    SS0(index) = Exporter_Prem(1,i) - Exporter_Prem_ROW(1,i); index = index+1;
    SS0(index) = TO(1,i) - TO_ROW(1,i); index = index+1;
end


end
