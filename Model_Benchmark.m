%%% Deterministic Model

%{

This file runs the Benchmark Model:
1. Calibrates a symmetric economy
2. Calibrates ROW
3. Calibrates SOE
4. Runs Dynare files

Two countries: SOE, ROW1
Multi Sector
Export dynamics: New exporters pay higher costs than continuing exporters
Factors of production: Labor, capital, materials & energy

6-Sector Model: AG, SS, Consumer Goods, Fuels, Ind. Supl, Capital and
Transp
%}


clear all;
clc;
close all;
%addpath c:\dynare\4.6.4\matlab;
addpath c:\dynare\4.5.7\matlab;

%% Parameters
global beta delta theta eta n sigma omega_c omega_k omega alpha theta_s fe fp f0 f1 xi_0 xi_1 tau S CO L ...
    phi_k kappa lambda sigma_m ...
    Exporter_SOE New_exporter_SOE EX_new_SOE Exporter_Prem_SOE TO_SOE ...
    Exporter_ROW New_exporter_ROW EX_new_ROW Exporter_Prem_ROW TO_ROW A TB_SOE

%% Definition of Parameters and initial values

S = 6; %number of sectors
CO = 2; %number of countries
L = [10; 10]; %Size
beta = 0.96;
eta   = 2.5; %Pareto distribution parameter
n     = [0.90; 0.90]; %Survival probability
delta = 0.10*ones(CO,S); %Capital depreciation
theta = 1.1; %Elasticity of substitution consumption goods
omega_c = 1/S*ones(CO,S); %Sector shares consumption

sigma = 1.1*ones(CO,S); %Elasticity of sustitution investment goods
omega_k = zeros(CO,S,S);
omega_k(1,:,:) = 1/S*ones(S,S); %Sector shares investment
omega_k(2,:,:) = 1/S*ones(S,S); 

%Shares CD production function
alpha = 0.31*ones(CO,S); %Capital share
kappa = 0.24*ones(CO,S); %Labor share

%Sector shares: Materials
lambda = zeros(CO,S,S);
lambda(1,:,:) = 1/S*ones(S,S);
lambda(2,:,:) = 1/S*ones(S,S);

sigma_m = 1.1*ones(CO,S); %Elasticity of substitution production of materials                
theta_s = 3.5*ones(CO,S); %Elasticity of substitution: Home-Foreign

%Home bias
omega = zeros(CO,S,CO);
omega(1,:,:) = 1/CO*ones(S,CO);           
omega(2,:,:) = 1/CO*ones(S,CO);
       
%Fixed Costs
fe = 5*ones(CO,1)/S/CO; %Entry cost
fp = ones(CO,S)/S/CO; %Fixed production cost
f0 = ones(CO,S)/S/CO; %Fixed export cost (new)
f1 = ones(CO,S)/S/CO; %Fixed export cost (old)

xi_0 = 1.15*ones(CO,S,CO); %From country CO to Country CO'
for m = 1:CO
xi_0(m,:,m) = ones(1,S);
end

xi_1 = 1.10*ones(CO,S,CO);  %Iceberg cost (continuation)
for m = 1:CO
xi_1(m,:,m) = ones(1,S);
end

tau = 1.40*ones(CO,S,CO); %Country that imposes the tariff first column
for m = 1:CO
tau(m,:,m) = ones(1,S);
end
          
phi_k = ones(CO,S); %Investment adjustment cost
phi_k(1,:) = ones(1,S); %Investment adjustment cost
phi_k(2,:) = ones(1,S); %Investment adjustment cost

%TFP
A = ones(CO,S);

%% Targets (Initial)

Exporter_ROW = 0.15*ones(1,S);
New_exporter_ROW = 0.50*ones(1,S);
EX_new_ROW = 0.15*ones(1,S);
Exporter_Prem_ROW = 1.5*ones(1,S);
TO_ROW = 0.20*ones(1,S);

%% Symmetric sectors: One Country

%Initial conditions
C = 10; W = 4; NT = 20; P = 1; K = 25; a_0 = 5; M = 3;
Xss0 = [C; W; NT; P; K; a_0; M];   
Xss1 = real(fsolve(@(x) Steady_State_Sym(x),Xss0));    

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new] = Steady_State_Sym(Xss1);


%% Calibrate Symmetric

while theta_s(1,1) > 3.0
theta_s(1,1) = theta_s(1,1) - 0.01; 
theta_s = theta_s(1,1)*ones(CO,S); %Elasticity of substitution: Home-Foreign

Xss0 = [Xss1; fp(1,1); f0(1,1); xi_0(1,1,2)];   
Xss1 = real(fsolve(@(x) Steady_State_SymCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new] = Steady_State_SymCalib(Xss1);

Xss0 = [C; W; NT; P(1,1); K(1,1); a_0(1,1); M(1,1)];
Xss1 = real(fsolve(@(x) Steady_State_Sym(x),Xss0));    

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T] = Steady_State_Sym(Xss1);

Xss1 = Xss0;
end

Xss1 = real(fsolve(@(x) Steady_State_Sym(x),Xss0));    

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T] = Steady_State_Sym(Xss1);

Xss1 = Xss0;

%% Calibrate Asymmetric Sectors
%% Step 1 
omega_c = [0.10, 0.25, 0.20, 0.15, 0.15, 0.15; ...
           0.10, 0.25, 0.20, 0.15, 0.15, 0.15]; 

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);

Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

%% Step 2
omega_c = [0.10, 0.35, 0.25, 0.075, 0.10, 0.075; ...
           0.10, 0.35, 0.25, 0.075, 0.10, 0.075]; 

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;
for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);


%% Step 3
omega_c = [0.11, 0.45, 0.33, 0.01, 0.07, 0.03; ...
           0.11, 0.45, 0.33, 0.01, 0.07, 0.03]; 

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;
for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

%% Step 4

lambda(1,:,:) = [0.38, 0.11, 0.20, 0.02, 0.24, 0.05; ...
                 0.06, 0.49, 0.11, 0.08, 0.22, 0.04; ...
                 0.47, 0.06, 0.34, 0.02, 0.10, 0.01; ...
                 0.33, 0.08, 0.07, 0.02, 0.48, 0.02;...
                 0.22, 0.12, 0.07, 0.04, 0.51, 0.04;...
                 0.03, 0.11, 0.05, 0.04, 0.38, 0.39];

                 
lambda(2,:,:) = [0.38, 0.11, 0.20, 0.02, 0.24, 0.05; ...
                 0.06, 0.49, 0.11, 0.08, 0.22, 0.04; ...
                 0.47, 0.06, 0.34, 0.02, 0.10, 0.01; ...
                 0.33, 0.08, 0.07, 0.02, 0.48, 0.02;...
                 0.22, 0.12, 0.07, 0.04, 0.51, 0.04;...
                 0.03, 0.11, 0.05, 0.04, 0.38, 0.39];

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);


%% Step 5

omega_k(1,:,:) = [0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15]; %Sector shares investment

omega_k(2,:,:) = [0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15;...
                  0.10, 0.30, 0.15, 0.15, 0.15, 0.15]; %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

%% Step 5

omega_k(1,:,:) = [0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25]; %Sector shares investment

omega_k(2,:,:) = [0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25;...
                  0.10, 0.35, 0.10, 0.10, 0.10, 0.25]; %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

%% Step 6

omega_k(1,:,:) = [0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25]; %Sector shares investment

omega_k(2,:,:) = [0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25;...
                  0.10, 0.50, 0.05, 0.05, 0.05, 0.25]; %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);


%% Step 7

omega_k(1,:,:) = [0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30]; %Sector shares investment

omega_k(2,:,:) = [0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30;...
                  0.05, 0.50, 0.05, 0.05, 0.05, 0.30]; %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

%% Step 8

omega_k(1,:,:) = [0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36]; %Sector shares investment

omega_k(2,:,:) = [0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36;...
                  0.05, 0.50, 0.03, 0.03, 0.03, 0.36]; %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);



%% Step 9

omega_k(1,:,:) = [0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419]; %Sector shares investment


omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment

Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = Xss1;

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);


%% Step 7

Exporter_ROW = [0.10 0.04 0.40 0.10 0.30 0.20];
New_exporter_ROW = [0.40 0.80 0.40 0.40 0.40 0.40];
EX_new_ROW = [0.15 0.40 0.15 0.15 0.15 0.15];
Exporter_Prem_ROW = [2.5 2.5 2.5 2.5  2.5 2.5];
TO_ROW = [0.10 0.04 0.20 0.20 0.40 0.40];

Exporter_SOE = [0.20 0.015 0.20 0.20 0.20 0.20];
New_exporter_SOE = [0.35 0.45 0.35 0.35 0.35 0.35];
EX_new_SOE = [0.075 0.15 0.075 0.075 0.075 0.075];
Exporter_Prem_SOE = [2.5 2.5 2.5 2.5 2.5 2.5];
TO_SOE = [0.35 0.025 0.15 0.55 0.40 0.85];
TB_SOE = [0.05 0.00 0.01 -0.01 -0.02 -0.02];

while Exporter_Prem_ROW(1,2) <= 3.75

Exporter_Prem_ROW = Exporter_Prem_ROW + 0.05;
Exporter_Prem_SOE = Exporter_Prem_ROW;

Xss0 = [C, W, NT];
index = 4;
for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSector(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSector(Xss1);


Xss0 = [C, W, NT];
index = 4;

for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymSectorCalib(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, P_M, X_M, M, MC_CD, CD, PI_CD_T, ...
    EX_1, EX_0, GO, GO_EX, GO_NoEX, EX_GO, TO, Exporter, New_exp, Exporter_Prem, EX_new, mu, D_M] = Steady_State_AsymSectorCalib(Xss1);

end

%% Asymmetric countries

%% Step 1

Xss0 = [P_C(1,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(1,1), W(1,1), NT(1,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(1,i); index = index + 1;
Xss0(index) = K(1,i); index = index + 1;
Xss0(index) = a_0(1,i); index = index + 1;
Xss0(index) = a_1(1,i); index = index + 1;
Xss0(index) = M(1,i); index = index + 1;
end
end

Xss1 = real(fsolve(@(x) Steady_State_Asym(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_Asym(Xss1);
Xss0 = Xss1;

%% Step 4: Tariffs

%tau(1,:,2) = [1.40, 1.40, 1.40, 1.30, 1.30, 1.25]; %Country that imposes the tariff first column

while tau(1,3,2) >= 1.30
    tau(1,3,2) = tau(1,3,2) - 0.001; 

    for i = 4:S
        tau(1,i,2) = tau(1,3,2);
    end

Xss0 = [P_C(2,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(2,1), W(2,1), NT(2,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(j,i); index = index + 1;
Xss0(index) = K(j,i); index = index + 1;
Xss0(index) = a_0(j,i); index = index + 1;
Xss0(index) = a_1(j,i); index = index + 1;
Xss0(index) = M(j,i); index = index + 1;
end
end

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

for i = 1:S-1
Xss0(index) = omega(1,i,1); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

end

%% Step 5: Tariffs

while tau(1,6,2) >= 1.25
    tau(1,6,2) = tau(1,6,2) - 0.001; 

Xss0 = [P_C(2,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(2,1), W(2,1), NT(2,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(j,i); index = index + 1;
Xss0(index) = K(j,i); index = index + 1;
Xss0(index) = a_0(j,i); index = index + 1;
Xss0(index) = a_1(j,i); index = index + 1;
Xss0(index) = M(j,i); index = index + 1;
end
end

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

for i = 1:S-1
Xss0(index) = omega(1,i,1); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

end


tau(1,:,2) = [1.40, 1.40, 1.40, 1.30, 1.30, 1.25];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;


tau(1,:,2) = [1.30, 1.40, 1.35, 1.10, 1.25, 1.25];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

while tau(1,2,2) > 1
    tau(1,2,2) = tau(1,2,2)-0.1;
Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;
end

tau(1,2,2) = 1;
Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

%% Trade Balance: Target Productivity

%% Step 1

Xss0 = [P_C(2,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(2,1), W(2,1), NT(2,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(j,i); index = index + 1;
Xss0(index) = K(j,i); index = index + 1;
Xss0(index) = a_0(j,i); index = index + 1;
Xss0(index) = a_1(j,i); index = index + 1;
Xss0(index) = M(j,i); index = index + 1;
end
end

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
end

for i = 1:S-1
Xss0(index) = omega(1,i,1); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

%% Step 2

TB_SOE = [0.10 -0.01 0.04 0.06 -0.12 -0.37];

while L(2,1) < 150

size = sum(GO(1,:))/(sum(GO(1,:))+sum(GO(2,:)));
L(2,1) = L(2,1) + 0.5;

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;
size = sum(GO(1,:))/(sum(GO(1,:))+sum(GO(2,:)));

end

%% Step 3
save Ini_ss0
load Ini_ss0

f0(2,:) = 1.5*f1(2,:);
xi_0(2,:,1) = 1.15*xi_1(2,:,1);

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

jj = 0.5;
omega(2,:,:) = jj*ones(S,CO);
while jj < 0.90
jj = jj +0.01;
    omega(2,:,:) = jj*ones(S,CO);

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

end


%% Data SOE

Exporter_SOE = [0.18 0.015 0.18 0.35 0.12 0.12];
Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

Exporter_SOE = [0.20 0.02 0.18 0.35 0.12 0.12];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

Exporter_Prem_SOE = [4.0 3.5 2.5 2.0 3.8 4.0];

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

Exporter_Prem_SOE = [3.5 3.5 2.5 2.0 3.8 4.0];
[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

Exporter_Prem_SOE = [3.0 3.5 2.5 2.0 3.8 4.0];
[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

TB_SOE = [0.10 -0.01 0.04 0.05 -0.15 -0.37];
Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

TO_SOE = [0.30 0.04 0.14 0.57 0.41 0.8];
TB_SOE = [0.20 -0.01 0.04 0.06 -0.20 -1.30];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;

while TB_SOE(1,1) < 0.25
TB_SOE(1,1) = TB_SOE(1,1) + 0.01;

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB(Xss1);
Xss0 = Xss1;
end

%% Calibrate New Exporters

New_exporter_SOE = [0.80 0.80 0.80 0.80 0.80 0.80];
EX_new_SOE = [0.80 0.80 0.80 0.80 0.80 0.80];

while New_exporter_SOE(1,3) > 0.35 

New_exporter_SOE(1,1) = New_exporter_SOE(1,1) - 0.01;
EX_new_SOE(1,1) = EX_new_SOE(1,1) - 0.015;

New_exporter_SOE(1,3) = New_exporter_SOE(1,3) - 0.01;
EX_new_SOE(1,3) = EX_new_SOE(1,3) - 0.015;

New_exporter_SOE(1,4) = New_exporter_SOE(1,4) - 0.01;
EX_new_SOE(1,4) = EX_new_SOE(1,4) - 0.015;

New_exporter_SOE(1,5) = New_exporter_SOE(1,5) - 0.01;
EX_new_SOE(1,5) = EX_new_SOE(1,5) - 0.015;

New_exporter_SOE(1,6) = New_exporter_SOE(1,6) - 0.01/2;
EX_new_SOE(1,6) = EX_new_SOE(1,6) - 0.015/2;

Xss0 = [P_C(2,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(2,1), W(2,1), NT(2,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(j,i); index = index + 1;
Xss0(index) = K(j,i); index = index + 1;
Xss0(index) = a_0(j,i); index = index + 1;
Xss0(index) = a_1(j,i); index = index + 1;
Xss0(index) = M(j,i); index = index + 1;
end
end

for i = 1:S
Xss0(index) = fp(1,i); index = index+1;
Xss0(index) = f0(1,i); index = index+1;
Xss0(index) = f1(1,i); index = index+1;
Xss0(index) = xi_0(1,i,2); index = index+1;
Xss0(index) = xi_1(1,i,2); index = index+1;
end

for i = 1:S-1
Xss0(index) = omega(1,i,1); index = index+1;
end

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

end

New_exporter_SOE = [0.30 0.38 0.32 0.20 0.30 0.35];
EX_new_SOE = [0.075 0.05 0.12 0.035 0.075 0.075];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;


%% 

size = sum(GO(1,:))/(sum(GO(1,:))+sum(GO(2,:)));

%% Check Moments

Moment = {'Share GO'; 'Share N'; 'TO'; 'Exporter'; 'New Exporter';'Exporter Premium';'EX_new/EX'; 'TB'};
Data_Sector_1 = [GO(1,1)/sum(GO(1,:)); N(1,1)/sum(N(1,:)); TO(1,1); Exporter(1,1); New_exp(1,1); Exporter_Prem(1,1); EX_0(1,1)/EX(1,1); TB(1,1)];
Data_Sector_2 = [GO(1,2)/sum(GO(1,:)); N(1,2)/sum(N(1,:)); TO(1,2); Exporter(1,2); New_exp(1,2); Exporter_Prem(1,2); EX_0(1,2)/EX(1,2); TB(1,2)];
Data_Sector_3 = [GO(1,3)/sum(GO(1,:)); N(1,3)/sum(N(1,:)); TO(1,3); Exporter(1,3); New_exp(1,3); Exporter_Prem(1,3); EX_0(1,3)/EX(1,3); TB(1,3)];
Data_Sector_4 = [GO(1,4)/sum(GO(1,:)); N(1,4)/sum(N(1,:)); TO(1,4); Exporter(1,4); New_exp(1,4); Exporter_Prem(1,4); EX_0(1,4)/EX(1,4); TB(1,4)];
Data_Sector_5 = [GO(1,5)/sum(GO(1,:)); N(1,5)/sum(N(1,:)); TO(1,5); Exporter(1,5); New_exp(1,5); Exporter_Prem(1,5); EX_0(1,5)/EX(1,5); TB(1,5)];
Data_Sector_6 = [GO(1,6)/sum(GO(1,:)); N(1,6)/sum(N(1,:)); TO(1,6); Exporter(1,6); New_exp(1,6); Exporter_Prem(1,6); EX_0(1,6)/EX(1,6); TB(1,6)];

T_SOE = table(Moment,Data_Sector_1, Data_Sector_2, Data_Sector_3, Data_Sector_4, Data_Sector_5, Data_Sector_6)

Moment = {'Share GO'; 'Share N'; 'TO'; 'Exporter'; 'New Exporter';'Exporter Premium';'EX_new/EX'; 'TB'; 'TradeShare'};
Data_Sector_1 = [GO(2,1)/sum(GO(2,:)); N(2,1)/sum(N(2,:)); TO(2,1); Exporter(2,1); New_exp(2,1); Exporter_Prem(2,1); EX_0(2,1)/EX(1,1); TB(2,1); (EX(2,1)+IM(2,1))/2/EX_T(2,1)];
Data_Sector_2 = [GO(2,2)/sum(GO(2,:)); N(2,2)/sum(N(2,:)); TO(2,2); Exporter(2,2); New_exp(2,2); Exporter_Prem(2,2); EX_0(2,2)/EX(1,2); TB(2,2); (EX(2,2)+IM(2,2))/2/EX_T(2,1)];
Data_Sector_3 = [GO(2,3)/sum(GO(2,:)); N(2,3)/sum(N(2,:)); TO(2,3); Exporter(2,3); New_exp(2,3); Exporter_Prem(2,3); EX_0(2,3)/EX(1,3); TB(2,3); (EX(2,3)+IM(2,3))/2/EX_T(2,1)];
Data_Sector_4 = [GO(2,4)/sum(GO(2,:)); N(2,4)/sum(N(2,:)); TO(2,4); Exporter(2,4); New_exp(2,4); Exporter_Prem(2,4); EX_0(2,4)/EX(1,4); TB(2,4); (EX(2,4)+IM(2,4))/2/EX_T(2,1)];
Data_Sector_5 = [GO(2,5)/sum(GO(2,:)); N(2,5)/sum(N(2,:)); TO(2,5); Exporter(2,5); New_exp(2,5); Exporter_Prem(2,5); EX_0(2,5)/EX(1,5); TB(2,5); (EX(2,5)+IM(2,5))/2/EX_T(2,1)];
Data_Sector_6 = [GO(2,6)/sum(GO(2,:)); N(2,6)/sum(N(2,:)); TO(2,6); Exporter(2,6); New_exp(2,6); Exporter_Prem(2,6); EX_0(2,6)/EX(1,6); TB(2,6); (EX(2,6)+IM(2,6))/2/EX_T(2,1)];

T_ROW = table(Moment,Data_Sector_1, Data_Sector_2, Data_Sector_3, Data_Sector_4, Data_Sector_5, Data_Sector_6)


size = sum(GO(1,:))/(sum(GO(1,:))+sum(GO(2,:)))

%% Macro Moments

Moment = {'TO/GDP'; 'C/GDP'; 'I/GDP'};
Data_SOE = [(EX_T(1,1)+IM_T(1,1))/(D_C(1,1)+sum(D_K(1,:))); D_C(1,1)/(D_C(1,1)+sum(D_K(1,:))); sum(D_K(1,:))/(D_C(1,1)+sum(D_K(1,:)))];
Data_ROW = [(EX_T(2,1)+IM_T(2,1))/(D_C(2,1)+sum(D_K(2,:))); D_C(2,1)/(D_C(2,1)+sum(D_K(2,:))); sum(D_K(2,:))/(D_C(2,1)+sum(D_K(2,:)))];
T_Macro = table(Moment,Data_SOE, Data_ROW)


save Ini_SS 

%% Dynare

%% Save Parameters
%Aggregate
Params_AG = [beta; eta; theta; S];

%Country
Params_CO = [L(1,1) L(2,1); n(1,1) n(2,1); ...
             fe(1,1) fe(2,1)];

        
%Sector-Country
Params_S = zeros(11,CO,S);

for m = 1:CO
    for j = 1:S
Params_S(:,m,j) = [sigma(m,j); ...
                   delta(m,j); ...
                   omega_c(m,j); ...
                   alpha(m,j); ...
                   theta_s(m,j); ...
                   fp(m,j); ...
                   f0(m,j); ...
                   f1(m,j); ...
                   phi_k(m,j); ...
                   kappa(m,j);...
                   sigma_m(m,j)];
    end
end

%Sector-Country-Investment
Params_I = omega_k;
Params_Omega = omega;
Params_lambda = lambda;


save Params_AG Params_AG
save Params_CO Params_CO
save Params_S Params_S
save Params_I Params_I
save Params_Omega Params_Omega
save Params_lambda Params_lambda
%% Save Variables

%Country
Xss0_CO = [P_C(1,1) P_C(2,1); Q Q; C(1,1) C(2,1); W(1,1) W(2,1); NT(1,1) NT(2,1); ...
           TauT(1,1) TauT(2,1); EX_T(1,1) EX_T(2,1); PI_T(1,1) PI_T(2,1); ...
           WL(1,1) WL(2,1); D_C(1,1) D_C(2,1); S*NE(1,1) S*NE(2,1);...
           IM_T(1,1) IM_T(2,1); Lp_T(1,1) Lp_T(2,1)];

%Country sector

Xss0_S =  zeros(44,CO,S);

for m = 1:CO
   for i = 1:S
    
       Xss0_S(:,m,i) = [P(m,i); K(m,i); P_K(m,i); R_K(m,i); X_C(m,i); I(m,i); D_K(m,i); ...
            X(m,i); D(m,i); MC(m,i); PI_0(m,i); a_p(m,i); a_0(m,i); a_1(m,i); ...
            Psi_p(m,i); Psi_0(m,i); Psi_1(m,i); n_p(m,i); n_0(m,i); n_1(m,i); EV_0(m,i); ...
            EV_1(m,i); NX(m,i); NX_1(m,i); NX_0(m,i); N_0(m,i); Lp(m,i); ...
            Lf(m,i); TauR(m,i); EX(m,i); PI(m,i); EV_inf(m,i); mu(m,i); ...
            f0(m,i); f1(m,i); fp(m,i); IM(m,i); P_M(m,i); mm(m,i); M(m,i); D_M(m,i);...
            MC_CD(m,i); CD(m,i); GO(m,i)];
       
   end
end

Xss0_K = X_K;
Xss0_M = X_M;
Tau0 = tau;
Ice0_0 = xi_0;
Ice1_0 = xi_1;

save Xss0_CO Xss0_CO
save Xss0_S Xss0_S
save Xss0_K Xss0_K
save Tau0 Tau0
save Ice0_0 Ice0_0
save Ice1_0 Ice1_0
save Xss0_M Xss0_M;


%% New Steady state
tau0 = tau;
jj = 1;
while tau(1,3,2)-1 > (tau0(1,3,2)-1)*0.50

jj = jj - 0.01;

tau(1,3,2) = (tau0(1,3,2)-1)*jj+1;
tau(1,1,2) = tau(1,3,2);

for i = 5:S
tau(1,i,2) = (tau0(1,i,2)-1)*jj+1;
end

Xss0 = [P_C(2,1),...
        C(1,1), W(1,1), NT(1,1),...
        C(2,1), W(2,1), NT(2,1)];
index = 8;

for j = 1:CO
for i = 1:S
Xss0(index) = P(j,i); index = index + 1;
Xss0(index) = K(j,i); index = index + 1;
Xss0(index) = a_0(j,i); index = index + 1;
Xss0(index) = a_1(j,i); index = index + 1;
Xss0(index) = M(j,i); index = index + 1;
end
end

Xss1 = real(fsolve(@(x) Steady_State_Asym(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_Asym(Xss1);
Xss0 = Xss1;

end

Xss1 = real(fsolve(@(x) Steady_State_Asym(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S] = Steady_State_Asym(Xss1);
Xss0 = Xss1;

for i = 5:S
tau(1,i,2) = tau(1,4,2);
end

tau(1,6,2) = 1.05;

Xss1 = real(fsolve(@(x) Steady_State_Asym(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S] = Steady_State_Asym(Xss1);
Xss0 = Xss1;

Moment = {'Share GO'; 'Share N'; 'TO'; 'Exporter'; 'New Exporter';'Exporter Premium';'EX_new/EX'; 'TB'};
Data_Sector_1 = [GO(1,1)/sum(GO(1,:)); N(1,1)/sum(N(1,:)); TO(1,1); Exporter(1,1); New_exp(1,1); Exporter_Prem(1,1); EX_0(1,1)/EX(1,1); TB(1,1)];
Data_Sector_2 = [GO(1,2)/sum(GO(1,:)); N(1,2)/sum(N(1,:)); TO(1,2); Exporter(1,2); New_exp(1,2); Exporter_Prem(1,2); EX_0(1,2)/EX(1,2); TB(1,2)];
Data_Sector_3 = [GO(1,3)/sum(GO(1,:)); N(1,3)/sum(N(1,:)); TO(1,3); Exporter(1,3); New_exp(1,3); Exporter_Prem(1,3); EX_0(1,3)/EX(1,3); TB(1,3)];
Data_Sector_4 = [GO(1,4)/sum(GO(1,:)); N(1,4)/sum(N(1,:)); TO(1,4); Exporter(1,4); New_exp(1,4); Exporter_Prem(1,4); EX_0(1,4)/EX(1,4); TB(1,4)];
Data_Sector_5 = [GO(1,5)/sum(GO(1,:)); N(1,5)/sum(N(1,:)); TO(1,5); Exporter(1,5); New_exp(1,5); Exporter_Prem(1,5); EX_0(1,5)/EX(1,5); TB(1,5)];
Data_Sector_6 = [GO(1,6)/sum(GO(1,:)); N(1,6)/sum(N(1,:)); TO(1,6); Exporter(1,6); New_exp(1,6); Exporter_Prem(1,6); EX_0(1,6)/EX(1,6); TB(1,6)];

T_SOE = table(Moment,round(Data_Sector_1,2), round(Data_Sector_2,2), round(Data_Sector_3,2), round(Data_Sector_4,2), round(Data_Sector_5,2), round(Data_Sector_6,2))


%Country
Xss1_CO = [P_C(1,1) P_C(2,1); Q Q; C(1,1) C(2,1); W(1,1) W(2,1); NT(1,1) NT(2,1); ...
           TauT(1,1) TauT(2,1); EX_T(1,1) EX_T(2,1); PI_T(1,1) PI_T(2,1); ...
           WL(1,1) WL(2,1); D_C(1,1) D_C(2,1); S*NE(1,1) S*NE(2,1);...
           IM_T(1,1) IM_T(2,1); Lp_T(1,1) Lp_T(2,1)];

%Country sector

Xss1_S =  zeros(44,CO,S);

for m = 1:CO
   for i = 1:S
    
       Xss1_S(:,m,i) = [P(m,i); K(m,i); P_K(m,i); R_K(m,i); X_C(m,i); I(m,i); D_K(m,i); ...
            X(m,i); D(m,i); MC(m,i); PI_0(m,i); a_p(m,i); a_0(m,i); a_1(m,i); ...
            Psi_p(m,i); Psi_0(m,i); Psi_1(m,i); n_p(m,i); n_0(m,i); n_1(m,i); EV_0(m,i); ...
            EV_1(m,i); NX(m,i); NX_1(m,i); NX_0(m,i); N_0(m,i); Lp(m,i); ...
            Lf(m,i); TauR(m,i); EX(m,i); PI(m,i); EV_inf(m,i); mu(m,i); ...
            f0(m,i); f1(m,i); fp(m,i); IM(m,i); P_M(m,i); mm(m,i); M(m,i); D_M(m,i);...
            MC_CD(m,i); CD(m,i); GO(m,i)];
       
   end
end

Xss1_K = X_K;
Xss1_M = X_M;
Tau1 = tau;
Ice0_1 = xi_0;
Ice1_1 = xi_1;

save Xss1_CO Xss1_CO
save Xss1_S Xss1_S
save Xss1_K Xss1_K
save Tau1 Tau1
save Ice0_1 Ice0_1
save Ice1_1 Ice1_1
save Xss1_M Xss1_M;

save Final_ss

%% Save Values

dynare Model_NS_Bond.mod
oo_Het = oo_;
M_Het = M_;
save oo_Het oo_Het
save M_Het M_Het


