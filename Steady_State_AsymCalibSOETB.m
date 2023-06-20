function [SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB] = Steady_State_AsymCalibSOETB(x)

global beta delta theta eta n sigma omega_c omega_k omega alpha theta_s fe fp f0 f1 xi_0 xi_1 tau S CO L  ...
     kappa sigma_m lambda Exporter_SOE New_exporter_SOE EX_new_SOE Exporter_Prem_SOE TO_SOE TB_SOE A
index = 1;

C = zeros(CO,1);
W = zeros(CO,1);
NT = zeros(CO,1);

P_C(2,1) = x(index); index = index+1;

for m=1:CO
C(m,1) = x(index); index = index+1;
W(m,1) = x(index); index = index+1;
NT(m,1) = x(index); index = index+1;
end

%Variables
P = zeros(CO,S);
K = zeros(CO,S);
P_K = zeros(CO,S);
R_K = zeros(CO,S);
X_C = zeros(CO,S); 
I  = zeros(CO,S);
D_K  = zeros(CO,S);
X_K = zeros(CO,S,S);
X = zeros(CO,S);
D = zeros(CO,S);
MC = zeros(CO,S);
PI_0 = zeros(CO,S);
a_p = zeros(CO,S);
a_0 = zeros(CO,S);
a_1 = zeros(CO,S);
Psi_p = zeros(CO,S);
Psi_0 = zeros(CO,S);
Psi_1 = zeros(CO,S);
n_p = zeros(CO,S);
n_0 = zeros(CO,S);
n_1 = zeros(CO,S);
EV_0 = zeros(CO,S); 
EV_1 = zeros(CO,S); 
N = zeros(CO,S); 
NE = zeros(CO,S); 
NX = zeros(CO,S); 
NX_1 = zeros(CO,S); 
NX_0 = zeros(CO,S); 
N_0 = zeros(CO,S); 
Lp = zeros(CO,S); 
Lf = zeros(CO,S); 
TauR = zeros(CO,S); 
EX = zeros(CO,S); 
EX_0 = zeros(CO,S); 
EX_1 = zeros(CO,S); 
PI = zeros(CO,S); 
EV_inf = zeros(CO,S);
mu = zeros(CO,S);
Exporter = zeros(CO,S);
New_exp  = zeros(CO,S);
IM  = zeros(CO,S);
TO  = zeros(CO,S);

P_M = zeros(CO,S);
mm = zeros(CO,S);
X_M = zeros(CO,S,S);
M = zeros(CO,S);
D_M = zeros(CO,S);
MC_CD = zeros(CO,S);

for j = 1:CO
for i = 1:S
    P(j,i) = x(index); index = index+1;
    K(j,i) = x(index); index = index+1;
    a_0(j,i) = x(index); index = index+1;
    a_1(j,i) = x(index); index = index+1;
    M(j,i) = x(index); index = index+1;
end
end

%Parameters

for i = 1:S
fp(1,i) = x(index); index = index+1;
f0(1,i) = x(index); index = index+1;
f1(1,i) = f0(1,i);
xi_0(1,i,2) = x(index); index = index+1;
xi_1(1,i,2) = xi_0(1,i,2);
end

for i = 1:S-1
omega(1,i,1)  = x(index); index = index+1;
end
omega(1,:,2) = 1-omega(1,:,1);

%Consumption price Normalization
P_C(1,1) = 1;
%Nominal expenditure in C
D_C = zeros(CO,1);
for m = 1:CO
D_C(m,1) = P_C(m,1)*C(m,1);
end

%Stochastic discount factor
Q = beta;

%Price of sector-materials
for m = 1:CO
for i = 1:S
for j = 1:S
    P_M(m,i) = P_M(m,i) + lambda(m,i,j)*P(m,j)^(1-sigma_m(m,i));
end
    P_M(m,i) = P_M(m,i)^(1/(1-sigma_m(m,i)));
end
end

%Price of sector-investment
for m = 1:CO
for i = 1:S
for j = 1:S
    P_K(m,i) = P_K(m,i) + omega_k(m,i,j)*P(m,j)^(1-sigma(m,i));
end
    P_K(m,i) = P_K(m,i)^(1/(1-sigma(m,i)));
end
end

for m=1:CO
for i=1:S
%Interest rate
R_K(m,i) = (1/beta - (1-delta(1,i)))*P_K(m,i);
%Demand for sector-intermediates: consumption
X_C(m,i) = omega_c(m,i)*P(m,i)^(-theta)*P_C(m,1)^(theta-1)*D_C(m,1);
%Investment
I(m,i) = K(m,i)-(1-delta(m,i))*K(m,i);
%Nominal expenditure on investment
D_K(m,i) = I(m,i)*P_K(m,i);
%Nominal expenditure on materials
D_M(m,i) = P_M(m,i)*M(m,i);
end
end

for m = 1:CO
for i = 1:S
for j = 1:S
    %Demand for sector-intermediates: sector-investment
    X_K(m,i,j) = omega_k(m,i,j)*P(m,j)^(-sigma(m,i))*P_K(m,i)^(sigma(m,i)-1)*D_K(m,i);
    %Demand for sector-intermediates: sector-materials
    X_M(m,i,j) = lambda(m,i,j)*P(m,j)^(-sigma_m(m,i))*P_M(m,i)^(sigma_m(m,i)-1)*D_M(m,i);
end
end
end


%Total demand for sector intermediates
for m = 1:CO
for i = 1:S
X(m,i) = X_C(m,i);
for j = 1:S
   X(m,i) = X(m,i) + X_K(m,j,i)+ X_M(m,j,i);
end
end
end


for m = 1:CO
for i = 1:S
%Nominal expenditure on sector intermediates
D(m,i) = P(m,i)*X(m,i);
%Heterogeneous firms
%Marginal cost
MC_CD(m,i) = 1/A(m,i)*(R_K(m,i)/alpha(m,i))^alpha(m,i)*(W(m,1)/kappa(m,i))^kappa(m,i)*(P_M(m,i)/(1-alpha(m,i)-kappa(m,i)))^(1-alpha(m,i)-kappa(m,i));
MC(m,i) = MC_CD(m,i);   
end
end

for m = 1:CO
for i = 1:S
%Profits constant
PI_0(m,i) = 1/theta_s(m,i)*(theta_s(m,i)*MC(m,i)/(theta_s(m,i)-1))^(1-theta_s(m,i));
%Marginal producer: productivity
a_p(m,i) = W(m,1)*(fp(m,i))/(PI_0(m,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m));
end
end

%Distributions
Psi_p = eta/(eta-1).*a_p.^(1-eta);
Psi_0 = eta/(eta-1).*a_0.^(1-eta);
Psi_1 = eta/(eta-1).*a_1.^(1-eta);

n_p = a_p.^(-eta);
n_0 = a_0.^(-eta);
n_1 = a_1.^(-eta);

for m = 1:CO
for i = 1:S
%Expected values
VV = [1,  -n(m,1)*Q*n_1(m,i), -n(m,1)*Q*(1-n_1(m,i));...
      0, 1-n(m,1)*Q*n_1(m,i), -n(m,1)*Q*(1-n_1(m,i));...
      -n(m,1)*Q*n_0(m,i), 0, 1-(1-n_0(m,i))*n(m,1)*Q];

if m == 1   
EE = [PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i)+xi_0(m,i,m+1)^(1-theta_s(m+1,i))*tau(m+1,i,m)^(-theta_s(m+1,i))*P(m+1,i)^(theta_s(m+1,i)-1)*D(m+1,i)*(omega(m+1,i,m))*Psi_0(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)-W(m,1)*f0(m,i)*n_0(m,i); 
      PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i)+xi_1(m,i,m+1)^(1-theta_s(m+1,i))*tau(m+1,i,m)^(-theta_s(m+1,i))*P(m+1,i)^(theta_s(m+1,i)-1)*D(m+1,i)*(omega(m+1,i,m))*Psi_1(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)-W(m,1)*f1(m,i)*n_1(m,i);
      PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)];

  
EVV = VV^(-1)*EE;
EV_0(m,i) = EVV(1,1);
EV_1(m,i) = EVV(2,1);
EV_inf(m,i) = EVV(3,1);

elseif m == 2
EE = [PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i)+xi_0(m,i,m-1)^(1-theta_s(m-1,i))*tau(m-1,i,m)^(-theta_s(m-1,i))*P(m-1,i)^(theta_s(m-1,i)-1)*D(m-1,i)*(omega(m-1,i,m))*Psi_0(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)-W(m,1)*f0(m,i)*n_0(m,i); 
      PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i)+xi_1(m,i,m-1)^(1-theta_s(m-1,i))*tau(m-1,i,m)^(-theta_s(m-1,i))*P(m-1,i)^(theta_s(m-1,i)-1)*D(m-1,i)*(omega(m-1,i,m))*Psi_1(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)-W(m,1)*f1(m,i)*n_1(m,i);
      PI_0(m,i)*(P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)*Psi_p(m,i))-W(m,1)*(fp(m,i))*n_p(m,i)];

EVV = VV^(-1)*EE;
EV_0(m,i) = EVV(1,1);
EV_1(m,i) = EVV(2,1);
EV_inf(m,i) = EVV(3,1); 
  
end
end

%Masses of firms

for m = 1:CO
for i = 1:S

N(m,i) = NT(m,1)/S;
NE(m,i) = (1-n(m,1))*N(m,i);

AA = [1, -1, -1, 0; ...
      -n(m,1)*n_1(m,i), 1, 0, 0; ...
      0, 0, 1, -n(m,1)*n_0(m,i); ...
      1, 0, 0, 1];

NN = [0; 0; n_0(m,i)*NE(m,i); N(m,i)];
%NN = [0; 0; 0; N(m,i)];
NXS = AA^(-1)*NN;

%Exporters
NX(m,i)   = NXS(1,1);
%Continuation exporters
NX_1(m,i) = NXS(2,1);
%New exporters
NX_0(m,i) = NXS(3,1);
%Non exporters
N_0(m,i)  = NXS(4,1);

%%Fraction of Exporters
Exporter(m,i) = 1 - N_0(m,i)/N(m,i);

%New Exporters
%New_exp(m,i) = (NX_0(m,i))./(Exporter(m,i).*N(m,i));
New_exp(m,i) = NX_0(m,i)/NX(m,i);
end
end

for m = 1:CO
for i = 1:S   
%Production labor
Lp(m,i) = (theta_s(m,i)-1)*PI_0(m,i)*(kappa(m,i))/W(m,1)*(N(m,i)*Psi_p(m,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m));
for j = 1:CO
    if j~=m
Lp(m,i) = Lp(m,i) + (theta_s(m,i)-1)*PI_0(m,i)*(kappa(m,i))/W(m,1)*(tau(j,i,m)^(-theta_s(j,i))*P(j,i)^(theta_s(j,i)-1)*D(j,i)*(omega(j,i,m))*(xi_1(m,i,j)^(1-theta_s(j,i))*NX_1(m,i)/n_1(m,i)*Psi_1(m,i)+ ...
           xi_0(m,i,j)^(1-theta_s(j,i))*NX_0(m,i)/n_0(m,i)*Psi_0(m,i)));
    end
end
end
end

for m = 1:CO
for i = 1:S
%Labor for fixed costs
Lf(m,i) = N(m,i)*n_p(m,i)*(fp(m,i)) + f0(m,i)*NX_0(m,i) + f1(m,i)*NX_1(m,i) + NE(m,i)*fe(m,1);
end
end


%Tariff revenue
TauR = zeros(CO,S);
%Exports
EX = zeros(CO,S);
EX_1 = zeros(CO,S);
EX_0 = zeros(CO,S);


for m = 1:CO
for i = 1:S
%Profits
PI(m,i) = PI_0(m,i)*(N(m,i)*Psi_p(m,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m)) - W(m,1)*Lf(m,i);
%Imports
IM(m,i) = 0;
for j = 1:CO
    if j~=m
TauR(m,i) = TauR(m,i) +(omega(m,i,j))*(tau(m,i,j)-1)*theta_s(m,i)*PI_0(j,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*tau(m,i,j)^(-theta_s(m,i))*...
            (xi_1(j,i,m)^(1-theta_s(m,i))*NX_1(j,i)/n_1(j,i)*Psi_1(j,i)+ ...
             xi_0(j,i,m)^(1-theta_s(m,i))*NX_0(j,i)/n_0(j,i)*Psi_0(j,i));
EX(m,i) = EX(m,i) + (omega(j,i,m))*theta_s(j,i)*PI_0(m,i)*P(j,i)^(theta_s(j,i)-1)*D(j,i)*tau(j,i,m)^(-theta_s(j,i))*...
            (xi_1(m,i,j)^(1-theta_s(j,i))*NX_1(m,i)/n_1(m,i)*Psi_1(m,i)+...
             xi_0(m,i,j)^(1-theta_s(j,i))*NX_0(m,i)/n_0(m,i)*Psi_0(m,i));
EX_1(m,i) = EX_1(m,i) + (omega(j,i,m))*theta_s(j,i)*PI_0(m,i)*P(j,i)^(theta_s(j,i)-1)*D(j,i)*tau(j,i,m)^(-theta_s(j,i))*...
            (xi_1(m,i,j)^(1-theta_s(j,i))*NX_1(m,i)/n_1(m,i)*Psi_1(m,i));
EX_0(m,i) = EX_0(m,i) + (omega(j,i,m))*theta_s(j,i)*PI_0(m,i)*P(j,i)^(theta_s(j,i)-1)*D(j,i)*tau(j,i,m)^(-theta_s(j,i))*...
            (xi_0(m,i,j)^(1-theta_s(j,i))*NX_0(m,i)/n_0(m,i)*Psi_0(m,i));        
PI(m,i) = PI(m,i) + PI_0(m,i)*(tau(j,i,m)^(-theta_s(j,i))*P(j,i)^(theta_s(j,i)-1)*D(j,i)*(omega(j,i,m))*...
          (xi_1(m,i,j)^(1-theta_s(m,i))*NX_1(m,i)/n_1(m,i)*Psi_1(m,i)+...
           xi_0(m,i,j)^(1-theta_s(m,i))*NX_0(m,i)/n_0(m,i)*Psi_0(m,i)));
IM(m,i) = IM(m,i) + (omega(m,i,j))*theta_s(m,i)*PI_0(j,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*tau(m,i,j)^(-theta_s(m,i))*...
            (xi_1(j,i,m)^(1-theta_s(m,i))*NX_1(j,i)/n_1(j,i)*Psi_1(j,i)+...
             xi_0(j,i,m)^(1-theta_s(m,i))*NX_0(j,i)/n_0(j,i)*Psi_0(j,i));       
    end
end

      
end
end

EX_S = zeros(CO,S,CO);
for m = 1:CO
for i = 1:S
    for j = 1:CO
    if j~=m
EX_S(m,i,j) = (omega(j,i,m))*theta_s(j,i)*PI_0(m,i)*P(j,i)^(theta_s(j,i)-1)*D(j,i)*tau(j,i,m)^(-theta_s(j,i))*...
            (xi_1(m,i,j)^(1-theta_s(j,i))*NX_1(m,i)/n_1(m,i)*Psi_1(m,i)+...
             xi_0(m,i,j)^(1-theta_s(j,i))*NX_0(m,i)/n_0(m,i)*Psi_0(m,i));
    end
    end
end
end


%Gross Output by sector
GO = zeros(CO,S);
GO_EX = zeros(CO,S);
GO_NoEX = zeros(CO,S);
EX_GO = zeros(CO,S);
TO = zeros(CO,S);
TB = zeros(CO,S);
for m = 1:CO
    for i = 1:S
        GO(m,i)    = theta_s(m,i)*PI_0(m,i)*N(m,i)*Psi_p(m,i)*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m) + EX(m,i);
        GO_EX(m,i) = theta_s(m,i)*PI_0(m,i)*(NX_1(m,i)*Psi_1(m,i)+NX_0(m,i)*Psi_0(m,i))*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m) + EX(m,i);
        GO_NoEX(m,i) = theta_s(m,i)*PI_0(m,i)*(N(m,i)*Psi_p(m,i)-NX_1(m,i)*Psi_1(m,i)-NX_0(m,i)*Psi_0(m,i))*P(m,i)^(theta_s(m,i)-1)*D(m,i)*omega(m,i,m);
        EX_GO(m,i) = EX(m,i)/GO(m,i);
        TO(m,i) = (EX(m,i)+IM(m,i))/GO(m,i);
        TB(m,i) = (EX(m,i)-IM(m,i))/GO(m,i);
        Exporter_Prem(m,i) = GO_EX(m,i)/GO(m,i)/Exporter(m,i); 
        EX_new(m,i) = EX_0(m,i)/EX(m,i);
    end
end

for m = 1:CO
for i = 1:S
mm(m,i) = (1-alpha(m,i)-kappa(m,i))/kappa(m,i)*W(m,1)/P_M(m,i)*Lp(m,i);
CD(m,i) = K(m,i)^alpha(m,i)*Lp(m,i)^kappa(m,i)*mm(m,i)^(1-alpha(m,i)-kappa(m,i));
PI_CD(m,i) = MC_CD(m,i)*CD(m,i) - R_K(m,i)*K(m,i) - W(m,1)*Lp(m,i) - P_M(m,i)*mm(m,i);
end
end

%%Totals

for m = 1:CO
TauT(m,1) = sum(TauR(m,:));
EX_T(m,1) = sum(EX(m,:));
PI_T(m,1) = sum(PI(m,:));
IM_T(m,1) = sum(IM(m,:));
Lp_T(m,1) = sum(Lp(m,:));
end

%% Check Walras: Consumer's Budget Constraint

for m = 1:CO
WL(m,1) = P_C(m,1)*C(m,1) - TauT(m,1) - PI_T(m,1) - W(m,1)*L(m,1);
    for i= 1:S
   WL(m,1) = WL(m,1) + P_K(m,i)*I(m,i) - R_K(m,i)*K(m,i);    
    end
end

for m=1:CO
for i= 1:S
   mu(m,i) = 1/C(m,1)/P_C(m,1)*P_K(m,i);    
end
end




%% Solve
index = 1;

for m = 1:CO
for i=1:S
%Marginal new exporter
SS0(index) = W(m,1)*f0(m,i)-n(m,1)*Q*(EV_1(m,i)-EV_inf(m,i)); 
for j = 1:CO
    if j~=m
   SS0(index) = SS0(index)  - PI_0(m,i)*a_0(m,i)*xi_0(m,i,j)^(1-theta_s(j,i))*D(j,i)*(omega(j,i,m))/(tau(j,i,m)^theta_s(j,i)*P(j,i)^(1-theta_s(j,i)));
    end
end
index = index + 1;

%Marginal old exporter
SS0(index) = W(m,1)*f1(m,i) - n(m,1)*Q*(EV_1(m,i)-EV_inf(m,i));
for j = 1:CO
    if j~=m
   SS0(index) = SS0(index)  - PI_0(m,i)*a_1(m,i)*xi_1(m,i,j)^(1-theta_s(j,i))*D(j,i)*(omega(j,i,m))/(tau(j,i,m)^theta_s(j,i)*P(j,i)^(1-theta_s(j,i)));
    end
end
index = index + 1;
   
end
end

for m = 1:CO
%Free entry
SS0(index) = W(m,1)*fe(m,1);
for i=1:S
    SS0(index) = SS0(index) - Q*EV_inf(m,i)/S;
end
index = index + 1;
end
%Sector prices

for m = 1:CO
for i = 1:S
    SS0(index) = P(m,i)^(1-theta_s(m,i)) - omega(m,i,m)*N(m,i)*(theta_s(m,i)*MC(m,i)/(theta_s(m,i)-1))^(1-theta_s(m,i))*Psi_p(m,i);
    for j = 1:CO
       if j~=m 
    
    SS0(index) = SS0(index) - ...
       (omega(m,i,j))*(tau(m,i,j)*theta_s(m,i)/(theta_s(m,i)-1)*MC(j,i))^(1-theta_s(m,i))*...
       (xi_1(j,i,m)^(1-theta_s(m,i))*NX_1(j,i)/n_1(j,i)*Psi_1(j,i)+ ...
        xi_0(j,i,m)^(1-theta_s(m,i))*NX_0(j,i)/n_0(j,i)*Psi_0(j,i)); 
       
       end   
    end
    index = index + 1;
end
end

for m = 1:CO
%Labor equilibrium
SS0(index) = L(m,1) - sum(Lp(m,:)) - sum(Lf(m,:)); index = index + 1;
end

for m = 1:CO
%Capital equilibrium
for i = 1:S
   SS0(index) = K(m,i) - alpha(m,i)/(kappa(m,i))*W(m,1)/R_K(m,i)*Lp(m,i); index = index+1; 
    end
end

for m = 1:CO
%Consumption prices
P_Ceq = 0;
for i = 1:S
   P_Ceq = P_Ceq + omega_c(m,i)*P(m,i)^(1-theta);
end
P_Ceq = P_Ceq^(1/(1-theta));
SS0(index) = P_C(m,1) - P_Ceq; index = index+1;
end

for m = 1:CO
for i = 1:S
   SS0(index) = M(m,i) - mm(m,i); index = index+1;
end
end

%%Targets
for i = 1:S
    SS0(index) = Exporter(1,i) - Exporter_SOE(1,i); index = index+1;
    SS0(index) = Exporter_Prem(1,i) - Exporter_Prem_SOE(1,i); index = index+1;
    SS0(index) = TO(1,i) - TO_SOE(1,i); index = index+1;
end

for i = 1:S-1
SS0(index) = TB(1,i) - TB_SOE(1,i); index = index+1;
end

%%

SS0(index) = EX_T(1,1) - IM_T(1,1); 

end

