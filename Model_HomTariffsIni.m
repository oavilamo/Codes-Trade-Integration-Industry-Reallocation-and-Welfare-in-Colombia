%% HomIniogeneizing tariffs

clear all
clc
close all

addpath c:\dynare\4.5.7\matlab;

load Ini_SS.mat

global beta delta theta eta n sigma omega_c omega_k omega alpha theta_s fe fp f0 f1 xi_0 xi_1 tau S CO L ...
    phi_k kappa lambda sigma_m ...
    Exporter_SOE New_exporter_SOE EX_new_SOE Exporter_Prem_SOE TO_SOE ...
    Exporter_ROW New_exporter_ROW EX_new_ROW Exporter_Prem_ROW TO_ROW A TB_SOE
%% New Steady state
tau0 = tau;

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

%% Tariff revenue
global TauT0
TauT0 = TauT(1,1)/IM_T(1,1);

Xss0 = [Xss0, tau(1,1,2)];
Xss1 = real(fsolve(@(x) Steady_State_AsymTau(x),Xss0));

[SS0, P_C, Q, C, W, NT, P, K, P_K, R_K, X_C, I, D_K, X_K, X, D, MC, PI_0, a_p, a_0, a_1, ...
    Psi_p, Psi_0, Psi_1, n_p, n_0, n_1, EV_0, EV_1, N, NE, NX, NX_1, NX_0, N_0, Lp, Lf, TauR, EX, ...
    PI, TauT, EX_T, PI_T, WL, D_C, EV_inf, mu, Exporter, New_exp, GO, EX_GO, IM, TO, IM_T, ...
    GO_EX, GO_NoEX, EX_0, EX_1, P_M, mm, X_M, M, D_M, MC_CD, CD, Lp_T, EX_S] = Steady_State_AsymTau(Xss1);
Xss0 = Xss1;


Moment = {'Share GO'; 'Share N'; 'TO'; 'Exporter'; 'New Exporter';'Exporter Premium';'EX_new/EX'; 'TB'};
Data_Sector_1 = [GO(1,1)/sum(GO(1,:)); N(1,1)/sum(N(1,:)); TO(1,1); Exporter(1,1); New_exp(1,1); Exporter_Prem(1,1); EX_0(1,1)/EX(1,1); TB(1,1)];
Data_Sector_2 = [GO(1,2)/sum(GO(1,:)); N(1,2)/sum(N(1,:)); TO(1,2); Exporter(1,2); New_exp(1,2); Exporter_Prem(1,2); EX_0(1,2)/EX(1,2); TB(1,2)];
Data_Sector_3 = [GO(1,3)/sum(GO(1,:)); N(1,3)/sum(N(1,:)); TO(1,3); Exporter(1,3); New_exp(1,3); Exporter_Prem(1,3); EX_0(1,3)/EX(1,3); TB(1,3)];
Data_Sector_4 = [GO(1,4)/sum(GO(1,:)); N(1,4)/sum(N(1,:)); TO(1,4); Exporter(1,4); New_exp(1,4); Exporter_Prem(1,4); EX_0(1,4)/EX(1,4); TB(1,4)];
Data_Sector_5 = [GO(1,5)/sum(GO(1,:)); N(1,5)/sum(N(1,:)); TO(1,5); Exporter(1,5); New_exp(1,5); Exporter_Prem(1,5); EX_0(1,5)/EX(1,5); TB(1,5)];
Data_Sector_6 = [GO(1,6)/sum(GO(1,:)); N(1,6)/sum(N(1,:)); TO(1,6); Exporter(1,6); New_exp(1,6); Exporter_Prem(1,6); EX_0(1,6)/EX(1,6); TB(1,6)];

T_SOE = table(Moment,Data_Sector_1, Data_Sector_2, Data_Sector_3, Data_Sector_4, Data_Sector_5, Data_Sector_6)


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

%% 
dynare Model_NSBondHom2.mod

oo_HetHomIni = oo_;
M_HetHomIni = M_;

save oo_HetHomIni oo_HetHomIni
save M_HetHomIni M_HetHomIni


