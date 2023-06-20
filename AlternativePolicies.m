%% Alternative Policies 
%% First run Model_Benchmark.m to generate initial steady state Ini_ss

load Ini_ss

%{
dum_Acc = 1 Accelerated
dum_Acc = 2 Delayed
dum_Acc = 3 Transitory
dum_Acc = 4 No Habit
dum_Acc = 5 Low Bond Adjustment Costs
dum_Acc = 6 No Bonds
%}

dum_Acc = 6; %Change the dummy for alternative policies

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

%% Dynare

if dum_Acc == 1
dynare Model_NS_BondHom.mod
oo_HetAcc = oo_;
M_HetAcc = M_;

save oo_HetAcc oo_HetAcc
save M_HetAcc M_HetAcc

elseif dum_Acc == 2
dynare Model_NS_BondDelay.mod
oo_HetDel = oo_;
M_HetDel = M_;

save oo_HetDel oo_HetDel
save M_HetDel M_HetDel

elseif dum_Acc == 3
dynare Model_NS_BondBack.mod
oo_HetBack = oo_;
M_HetBack = M_;

save oo_HetBack oo_HetBack
save M_HetBack M_HetBack

elseif dum_Acc == 4
dynare Model_NS_BondNoHab.mod
oo_HetNoHab = oo_;
M_HetNoHab = M_;

save oo_HetNoHab oo_HetNoHab
save M_HetNoHab M_HetNoHab

elseif dum_Acc == 5
dynare Model_NS_BondNoBAC.mod
oo_HetNoBAC = oo_;
M_HetNoBAC = M_;

save oo_HetNoBAC oo_HetNoBAC
save M_HetNoBAC M_HetNoBAC

elseif dum_Acc == 6

dynare Model_NS_NoBond.mod
oo_HetNoBond = oo_;
M_HetNoBond = M_;

save oo_HetNoBond oo_HetNoBond
save M_HetNoBond M_HetNoBond

end
