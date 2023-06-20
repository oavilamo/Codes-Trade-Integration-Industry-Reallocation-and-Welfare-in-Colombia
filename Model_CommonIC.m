%%% Deterministic Model

%{
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

%%Dummies for policies

%% Parameters

load Ini_ss.mat

%% 

global beta delta theta eta n sigma omega_c omega_k omega alpha theta_s fe fp f0 f1 xi_0 xi_1 tau S CO L ...
    phi_k kappa lambda sigma_m ...
    Exporter_SOE New_exporter_SOE EX_new_SOE Exporter_Prem_SOE TO_SOE ...
    Exporter_ROW New_exporter_ROW EX_new_ROW Exporter_Prem_ROW TO_ROW A TB_SOE

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

%% Original omega_c, omega_k

omega_c = [0.11, 0.45, 0.33, 0.01, 0.07, 0.03; ...
           0.11, 0.45, 0.33, 0.01, 0.07, 0.03]; 

omega_k(1,:,:) = [0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419;...
                  0.04, 0.52, 0.001/2, 0.001/2, 0.02, 0.419]; %Sector shares investment


omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

%% 

omega_c = [0.1 0.43 0.30 0.01 0.06 0.10;
           0.1 0.43 0.30 0.01 0.06 0.10];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;


omega_c = [0.09 0.46 0.25 0.01 0.06 0.13;
           0.09 0.46 0.25 0.01 0.06 0.13];

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

%% 
omega_k(1,:,:) = [0.09 0.46 0.10 0.01 0.06 0.28;
                  0.09 0.46 0.10 0.01 0.06 0.28;
                  0.09 0.46 0.10 0.01 0.06 0.28;
                  0.09 0.46 0.10 0.01 0.06 0.28;
                  0.09 0.46 0.10 0.01 0.06 0.28;
                  0.09 0.46 0.10 0.01 0.06 0.28]; %Sector shares investment


omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment


Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;


%% 
omega_k(1,:,:) = [0.09 0.46 0.15 0.01 0.06 0.23;
                  0.09 0.46 0.15 0.01 0.06 0.23;
                  0.09 0.46 0.15 0.01 0.06 0.23;
                  0.09 0.46 0.15 0.01 0.06 0.23;
                  0.09 0.46 0.15 0.01 0.06 0.23;
                  0.09 0.46 0.15 0.01 0.06 0.23]; %Sector shares investment


omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment


Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

%% 

while omega_k(1,1,3) < 0.25
    omega_k(1,:,3) = omega_k(1,:,3) + 0.01/2;
    omega_k(1,:,6) = omega_k(1,:,6) - 0.01/2;

omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment


Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;
    
end




%%

omega_k(1,:,:) = [0.09 0.46 0.25 0.01 0.06 0.13;
                  0.09 0.46 0.25 0.01 0.06 0.13;
                  0.09 0.46 0.25 0.01 0.06 0.13;
                  0.09 0.46 0.25 0.01 0.06 0.13;
                  0.09 0.46 0.25 0.01 0.06 0.13;
                  0.09 0.46 0.25 0.01 0.06 0.13]; %Sector shares investment


omega_k(2,:,:) = omega_k(1,:,:); %Sector shares investment


Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

%% 


Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

%% 

size = sum(GO(1,:))/(sum(GO(1,:))+sum(GO(2,:)));

%% Update trade

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

TO_SOE = 0.925*TO_SOE;
TB_SOE = 0.925*TB_SOE;

Xss1 = real(fsolve(@(x) Steady_State_AsymCalibSOETB2(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S, Exporter_Prem, TB]= Steady_State_AsymCalibSOETB2(Xss1);
Xss0 = Xss1;

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


save Ini_SSLE 

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

save Final_ssLE

%% 

dynare Model_NS_Bond.mod
oo_HetIC = oo_;
M_HetIC = M_;
save oo_HetIC oo_HetIC
save M_HetIC M_HetIC

