%% Welfare

load oo_Het.mat
load M_Het.mat
load oo_HetAcc.mat
load M_HetAcc.mat
load oo_HetDel.mat
load M_HetDel.mat
load oo_HetBack.mat
load M_HetBack.mat
load oo_HetNoHab.mat
load M_HetNoHab.mat
load oo_HetNoBAC.mat
load M_HetNoBAC.mat
load oo_HetNoBond.mat
load M_HetNoBond.mat
load oo_HetNOIO
load M_HetNOIO 
load oo_HetIC
load M_HetIC 
load oo_HetStatic.mat
load M_HetStatic.mat
load oo_HetStaticFA.mat
load M_HetStaticFA.mat

load oo_TwoS
load M_TwoS

load oo_sectorAG.mat
load oo_sectorMM1.mat
load oo_sectorMM3.mat
load oo_sectorMM4.mat
load M_sectorAG.mat
load M_sectorMM1.mat
load M_sectorMM3.mat
load M_sectorMM4.mat
load oo_HetHomIni.mat
load M_HetHomIni.mat
load oo_Hom.mat
load M_Hom.mat


%% Consumption
C_SOE_Het     =  oo_Het.endo_simul(strmatch('C_SOE',M_Het.endo_names,'exact'),1:end)';
C_SOE_HetAcc  =  oo_HetAcc.endo_simul(strmatch('C_SOE',M_HetAcc.endo_names,'exact'),1:end)';
C_SOE_HetDel  =  oo_HetDel.endo_simul(strmatch('C_SOE',M_HetDel.endo_names,'exact'),1:end)';
C_SOE_HetBack =  oo_HetBack.endo_simul(strmatch('C_SOE',M_HetBack.endo_names,'exact'),1:end)';
C_SOE_HetNoHab =  oo_HetNoHab.endo_simul(strmatch('C_SOE',M_HetNoHab.endo_names,'exact'),1:end)';
C_SOE_HetNoBAC =  oo_HetNoBAC.endo_simul(strmatch('C_SOE',M_HetNoBAC.endo_names,'exact'),1:end)';
C_SOE_HetNoBond =  oo_HetNoBond.endo_simul(strmatch('C_SOE',M_HetNoBond.endo_names,'exact'),1:end)';
C_SOE_HetNOIO  =  oo_HetNOIO.endo_simul(strmatch('C_SOE',M_HetNOIO.endo_names,'exact'),1:end)';
C_SOE_HetIC  =  oo_HetIC.endo_simul(strmatch('C_SOE',M_HetIC.endo_names,'exact'),1:end)';
C_SOE_HetStatic  =  oo_HetStatic.endo_simul(strmatch('C_SOE',M_HetStatic.endo_names,'exact'),1:end)';
C_SOE_HetStaticFA  =  oo_HetStaticFA.endo_simul(strmatch('C_SOE',M_HetStaticFA.endo_names,'exact'),1:end)';
C_SOE_TwoS  =  oo_TwoS.endo_simul(strmatch('C_SOE',M_TwoS.endo_names,'exact'),1:end)';
C_SOE_sectorAG  =  oo_sectorAG.endo_simul(strmatch('C_SOE',M_sectorAG.endo_names,'exact'),1:end)';
C_SOE_sectorMM1  =  oo_sectorMM1.endo_simul(strmatch('C_SOE',M_sectorMM1.endo_names,'exact'),1:end)';
C_SOE_sectorMM3 =  oo_sectorMM3.endo_simul(strmatch('C_SOE',M_sectorMM3.endo_names,'exact'),1:end-1)';
C_SOE_sectorMM4 =  oo_sectorMM4.endo_simul(strmatch('C_SOE',M_sectorMM4.endo_names,'exact'),1:end-1)';
C_SOE_HetHomIni  =  oo_HetHomIni.endo_simul(strmatch('C_SOE',M_HetHomIni.endo_names,'exact'),1:end)';
C_SOE_Hom  =  oo_Hom.endo_simul(strmatch('C_SOE',M_Hom.endo_names,'exact'),1:end)';


%Utility: U = (C{t}-phi_c*C{t-1})^(1-sigma_c)/(1-sigma_c)

sigma_c = 2;
phi_c = 0.75;

U_SOE_Het(1,1)    = (C_SOE_Het(1,1)-phi_c*C_SOE_Het(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetAcc(1,1) = (C_SOE_HetAcc(1,1)-phi_c*C_SOE_HetAcc(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetDel(1,1) = (C_SOE_HetDel(1,1)-phi_c*C_SOE_HetDel(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetBack(1,1) = (C_SOE_HetBack(1,1)-phi_c*C_SOE_HetBack(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoHab(1,1) = (C_SOE_HetNoHab(1,1)-0*C_SOE_HetNoHab(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoBAC(1,1) = (C_SOE_HetNoBAC(1,1)-phi_c*C_SOE_HetNoBAC(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoBond(1,1) = (C_SOE_HetNoBond(1,1)-phi_c*C_SOE_HetNoBond(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNOIO(1,1) = (C_SOE_HetNOIO(1,1)-phi_c*C_SOE_HetNOIO(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetIC(1,1) = (C_SOE_HetIC(1,1)-phi_c*C_SOE_HetIC(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetStatic(1,1) = (C_SOE_HetStatic(1,1)-phi_c*C_SOE_HetStatic(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetStaticFA(1,1) = (C_SOE_HetStaticFA(1,1)-phi_c*C_SOE_HetStaticFA(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_TwoS(1,1) = (C_SOE_TwoS(1,1)-phi_c*C_SOE_TwoS(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorAG(1,1) = (C_SOE_sectorAG(1,1)-phi_c*C_SOE_sectorAG(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM1(1,1) = (C_SOE_sectorMM1(1,1)-phi_c*C_SOE_sectorMM1(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM3(1,1) = (C_SOE_sectorMM3(1,1)-phi_c*C_SOE_sectorMM3(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM4(1,1) = (C_SOE_sectorMM4(1,1)-phi_c*C_SOE_sectorMM4(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetHomIni(1,1) = (C_SOE_HetHomIni(1,1)-phi_c*C_SOE_HetHomIni(1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_Hom(1,1) = (C_SOE_Hom(1,1)-phi_c*C_SOE_Hom(1,1))^(1-sigma_c)/(1-sigma_c);

for i = 2:500
U_SOE_Het(i,1)    = (C_SOE_Het(i,1)-phi_c*C_SOE_Het(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetAcc(i,1) = (C_SOE_HetAcc(i,1)-phi_c*C_SOE_HetAcc(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetDel(i,1) = (C_SOE_HetDel(i,1)-phi_c*C_SOE_HetDel(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetBack(i,1) = (C_SOE_HetBack(i,1)-phi_c*C_SOE_HetBack(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoHab(i,1) = (C_SOE_HetNoHab(i,1)-0*C_SOE_HetNoHab(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoBAC(i,1) = (C_SOE_HetNoBAC(i,1)-phi_c*C_SOE_HetNoBAC(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNoBond(i,1) = (C_SOE_HetNoBond(i,1)-phi_c*C_SOE_HetNoBond(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetNOIO(i,1) = (C_SOE_HetNOIO(i,1)-phi_c*C_SOE_HetNOIO(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetIC(i,1) = (C_SOE_HetIC(i,1)-phi_c*C_SOE_HetIC(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetStatic(i,1) = (C_SOE_HetStatic(i,1)-phi_c*C_SOE_HetStatic(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetStaticFA(i,1) = (C_SOE_HetStaticFA(i,1)-phi_c*C_SOE_HetStaticFA(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_TwoS(i,1) = (C_SOE_TwoS(i,1)-phi_c*C_SOE_TwoS(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorAG(i,1) = (C_SOE_sectorAG(i,1)-phi_c*C_SOE_sectorAG(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM1(i,1) = (C_SOE_sectorMM1(i,1)-phi_c*C_SOE_sectorMM1(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM3(i,1) = (C_SOE_sectorMM3(i,1)-phi_c*C_SOE_sectorMM3(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_sectorMM4(i,1) = (C_SOE_sectorMM4(i,1)-phi_c*C_SOE_sectorMM4(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_HetHomIni(i,1) = (C_SOE_HetHomIni(i,1)-phi_c*C_SOE_HetHomIni(i-1,1))^(1-sigma_c)/(1-sigma_c);
U_SOE_Hom(i,1) = (C_SOE_Hom(i,1)-phi_c*C_SOE_Hom(i-1,1))^(1-sigma_c)/(1-sigma_c);

end

beta = M_HetAcc.params(strmatch('beta',M_HetAcc.param_names,'exact'),1);

for i = 1:500
   Wel_SOE_Het(i,1) =  beta^(i-1)*U_SOE_Het(i,1);
   Wel_SOE_ini_Het(i,1) = beta^(i-1)*U_SOE_Het(1,1);
   Wel_SOE_HetAcc(i,1) =  beta^(i-1)*U_SOE_HetAcc(i,1);
   Wel_SOE_ini_HetAcc(i,1) = beta^(i-1)*U_SOE_HetAcc(1,1);
   Wel_SOE_HetDel(i,1) =  beta^(i-1)*U_SOE_HetDel(i,1);
   Wel_SOE_ini_HetDel(i,1) = beta^(i-1)*U_SOE_HetDel(1,1);   
   Wel_SOE_HetBack(i,1) =  beta^(i-1)*U_SOE_HetBack(i,1);
   Wel_SOE_ini_HetBack(i,1) = beta^(i-1)*U_SOE_HetBack(1,1);   
   Wel_SOE_HetNoHab(i,1) =  beta^(i-1)*U_SOE_HetNoHab(i,1);
   Wel_SOE_ini_HetNoHab(i,1) = beta^(i-1)*U_SOE_HetNoHab(1,1);
   Wel_SOE_HetNoBAC(i,1) =  beta^(i-1)*U_SOE_HetNoBAC(i,1);
   Wel_SOE_ini_HetNoBAC(i,1) = beta^(i-1)*U_SOE_HetNoBAC(1,1);
   Wel_SOE_HetNoBond(i,1) =  beta^(i-1)*U_SOE_HetNoBond(i,1);
   Wel_SOE_ini_HetNoBond(i,1) = beta^(i-1)*U_SOE_HetNoBond(1,1);
   Wel_SOE_HetNOIO(i,1) =  beta^(i-1)*U_SOE_HetNOIO(i,1);
   Wel_SOE_ini_HetNOIO(i,1) = beta^(i-1)*U_SOE_HetNOIO(1,1);
   Wel_SOE_HetIC(i,1) =  beta^(i-1)*U_SOE_HetIC(i,1);
   Wel_SOE_ini_HetIC(i,1) = beta^(i-1)*U_SOE_HetIC(1,1);
   Wel_SOE_HetStatic(i,1) =  beta^(i-1)*U_SOE_HetStatic(i,1);
   Wel_SOE_ini_HetStatic(i,1) = beta^(i-1)*U_SOE_HetStatic(1,1);
   Wel_SOE_HetStaticFA(i,1) =  beta^(i-1)*U_SOE_HetStaticFA(i,1);
   Wel_SOE_ini_HetStaticFA(i,1) = beta^(i-1)*U_SOE_HetStaticFA(1,1);
   Wel_SOE_TwoS(i,1) =  beta^(i-1)*U_SOE_TwoS(i,1);
   Wel_SOE_ini_TwoS(i,1) = beta^(i-1)*U_SOE_TwoS(1,1);
   Wel_SOE_Het(i,1) =  beta^(i-1)*U_SOE_Het(i,1);
   Wel_SOE_ini_Het(i,1) = beta^(i-1)*U_SOE_Het(1,1);
   Wel_SOE_sectorAG(i,1) =  beta^(i-1)*U_SOE_sectorAG(i,1);
   Wel_SOE_ini_sectorAG(i,1) = beta^(i-1)*U_SOE_sectorAG(1,1);
   Wel_SOE_sectorMM1(i,1) =  beta^(i-1)*U_SOE_sectorMM1(i,1);
   Wel_SOE_ini_sectorMM1(i,1) = beta^(i-1)*U_SOE_sectorMM1(1,1);   
   Wel_SOE_sectorMM3(i,1) =  beta^(i-1)*U_SOE_sectorMM3(i,1);
   Wel_SOE_ini_sectorMM3(i,1) = beta^(i-1)*U_SOE_sectorMM3(1,1);   
   Wel_SOE_sectorMM4(i,1) =  beta^(i-1)*U_SOE_sectorMM4(i,1);
   Wel_SOE_ini_sectorMM4(i,1) = beta^(i-1)*U_SOE_sectorMM4(1,1);   
   Wel_SOE_HetHomIni(i,1) =  beta^(i-1)*U_SOE_HetHomIni(i,1);
   Wel_SOE_ini_HetHomIni(i,1) = beta^(i-1)*U_SOE_HetHomIni(1,1);
   Wel_SOE_Hom(i,1) =  beta^(i-1)*U_SOE_Hom(i,1);
   Wel_SOE_ini_Hom(i,1) = beta^(i-1)*U_SOE_Hom(1,1);

end

PV_Wel_SOE_Het     = sum(Wel_SOE_Het);
PV_Wel_SOE_ini_Het = sum(Wel_SOE_ini_Het);
Delta_Wel_SOE_Het  = (PV_Wel_SOE_Het/PV_Wel_SOE_ini_Het)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_Het_ss = log(C_SOE_Het(end,1))-log(C_SOE_Het(1,1));

PV_Wel_SOE_HetAcc     = sum(Wel_SOE_HetAcc);
PV_Wel_SOE_ini_HetAcc = sum(Wel_SOE_ini_HetAcc);
Delta_Wel_SOE_HetAcc  = (PV_Wel_SOE_HetAcc/PV_Wel_SOE_ini_HetAcc)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetAcc_ss = log(C_SOE_HetAcc(end,1))-log(C_SOE_HetAcc(1,1));

PV_Wel_SOE_HetDel     = sum(Wel_SOE_HetDel);
PV_Wel_SOE_ini_HetDel= sum(Wel_SOE_ini_HetDel);
Delta_Wel_SOE_HetDel  = (PV_Wel_SOE_HetDel/PV_Wel_SOE_ini_HetDel)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetDel_ss = log(C_SOE_HetDel(end,1))-log(C_SOE_HetDel(1,1));

PV_Wel_SOE_HetBack     = sum(Wel_SOE_HetBack);
PV_Wel_SOE_ini_HetBack= sum(Wel_SOE_ini_HetBack);
Delta_Wel_SOE_HetBack  = (PV_Wel_SOE_HetBack/PV_Wel_SOE_ini_HetBack)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetBack_ss = log(C_SOE_HetBack(end,1))-log(C_SOE_HetBack(1,1));

PV_Wel_SOE_HetNoHab     = sum(Wel_SOE_HetNoHab);
PV_Wel_SOE_ini_HetNoHab= sum(Wel_SOE_ini_HetNoHab);
Delta_Wel_SOE_HetNoHab  = (PV_Wel_SOE_HetNoHab/PV_Wel_SOE_ini_HetNoHab)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetNoHab_ss = log(C_SOE_HetNoHab(end,1))-log(C_SOE_HetNoHab(1,1));

PV_Wel_SOE_HetNoBAC     = sum(Wel_SOE_HetNoBAC);
PV_Wel_SOE_ini_HetNoBAC= sum(Wel_SOE_ini_HetNoBAC);
Delta_Wel_SOE_HetNoBAC  = (PV_Wel_SOE_HetNoBAC/PV_Wel_SOE_ini_HetNoBAC)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetNoBAC_ss = log(C_SOE_HetNoBAC(end,1))-log(C_SOE_HetNoBAC(1,1));

PV_Wel_SOE_HetNoBond     = sum(Wel_SOE_HetNoBond);
PV_Wel_SOE_ini_HetNoBond = sum(Wel_SOE_ini_HetNoBond);
Delta_Wel_SOE_HetNoBond  = (PV_Wel_SOE_HetNoBond/PV_Wel_SOE_ini_HetNoBond)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetNoBond_ss = log(C_SOE_HetNoBond(end,1))-log(C_SOE_HetNoBond(1,1));

PV_Wel_SOE_HetNoIO     = sum(Wel_SOE_HetNOIO);
PV_Wel_SOE_ini_HetNoIO = sum(Wel_SOE_ini_HetNOIO);
Delta_Wel_SOE_HetNoIO  = (PV_Wel_SOE_HetNoIO/PV_Wel_SOE_ini_HetNoIO)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetNoIO_ss = log(C_SOE_HetNOIO(end,1))-log(C_SOE_HetNOIO(1,1));

PV_Wel_SOE_HetIC     = sum(Wel_SOE_HetIC);
PV_Wel_SOE_ini_HetIC = sum(Wel_SOE_ini_HetIC);
Delta_Wel_SOE_HetIC  = (PV_Wel_SOE_HetIC/PV_Wel_SOE_ini_HetIC)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetIC_ss = log(C_SOE_HetIC(end,1))-log(C_SOE_HetIC(1,1));

PV_Wel_SOE_HetStatic     = sum(Wel_SOE_HetStatic);
PV_Wel_SOE_ini_HetStatic = sum(Wel_SOE_ini_HetStatic);
Delta_Wel_SOE_HetStatic  = (PV_Wel_SOE_HetStatic/PV_Wel_SOE_ini_HetStatic)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetStatic_ss = log(C_SOE_HetStatic(end,1))-log(C_SOE_HetStatic(1,1));

PV_Wel_SOE_HetStaticFA     = sum(Wel_SOE_HetStaticFA);
PV_Wel_SOE_ini_HetStaticFA = sum(Wel_SOE_ini_HetStaticFA);
Delta_Wel_SOE_HetStaticFA  = (PV_Wel_SOE_HetStaticFA/PV_Wel_SOE_ini_HetStaticFA)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetStaticFA_ss = log(C_SOE_HetStaticFA(end,1))-log(C_SOE_HetStaticFA(1,1));

PV_Wel_SOE_TwoS     = sum(Wel_SOE_TwoS);
PV_Wel_SOE_ini_TwoS = sum(Wel_SOE_ini_TwoS);
Delta_Wel_SOE_TwoS  = (PV_Wel_SOE_TwoS/PV_Wel_SOE_ini_TwoS)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_TwoS_ss = log(C_SOE_TwoS(end,1))-log(C_SOE_TwoS(1,1));

PV_Wel_SOE_sectorAG     = sum(Wel_SOE_sectorAG);
PV_Wel_SOE_ini_sectorAG = sum(Wel_SOE_ini_sectorAG);
Delta_Wel_SOE_sectorAG  = (PV_Wel_SOE_sectorAG/PV_Wel_SOE_ini_sectorAG)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_sectorAG_ss = log(C_SOE_sectorAG(end,1))-log(C_SOE_sectorAG(1,1));

PV_Wel_SOE_sectorMM1     = sum(Wel_SOE_sectorMM1);
PV_Wel_SOE_ini_sectorMM1= sum(Wel_SOE_ini_sectorMM1);
Delta_Wel_SOE_sectorMM1  = (PV_Wel_SOE_sectorMM1/PV_Wel_SOE_ini_sectorMM1)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_sectorMM1_ss = log(C_SOE_sectorMM1(end,1))-log(C_SOE_sectorMM1(1,1));

PV_Wel_SOE_sectorMM3     = sum(Wel_SOE_sectorMM3);
PV_Wel_SOE_ini_sectorMM3= sum(Wel_SOE_ini_sectorMM3);
Delta_Wel_SOE_sectorMM3  = (PV_Wel_SOE_sectorMM3/PV_Wel_SOE_ini_sectorMM3)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_sectorMM3_ss = log(C_SOE_sectorMM3(end,1))-log(C_SOE_sectorMM3(1,1));

PV_Wel_SOE_sectorMM4     = sum(Wel_SOE_sectorMM4);
PV_Wel_SOE_ini_sectorMM4= sum(Wel_SOE_ini_sectorMM4);
Delta_Wel_SOE_sectorMM4  = (PV_Wel_SOE_sectorMM4/PV_Wel_SOE_ini_sectorMM4)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_sectorMM4_ss = log(C_SOE_sectorMM4(end,1))-log(C_SOE_sectorMM4(1,1));

PV_Wel_SOE_HetHomIni     = sum(Wel_SOE_HetHomIni);
PV_Wel_SOE_ini_HetHomIni = sum(Wel_SOE_ini_HetHomIni);
Delta_Wel_SOE_HetHomIni  = (PV_Wel_SOE_HetHomIni/PV_Wel_SOE_ini_HetHomIni)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_HetHomIni_ss = log(C_SOE_HetHomIni(end,1))-log(C_SOE_HetHomIni(1,1));

PV_Wel_SOE_Hom     = sum(Wel_SOE_Hom);
PV_Wel_SOE_ini_Hom = sum(Wel_SOE_ini_Hom);
Delta_Wel_SOE_Hom  = (PV_Wel_SOE_Hom/PV_Wel_SOE_ini_Hom)^(1/(1-sigma_c))-1;
Delta_Wel_SOE_Hom_ss = log(C_SOE_Hom(end,1))-log(C_SOE_Hom(1,1));

Moment = {'Benchmark'; 'AG'; 'MM Cons'; 'MM Ind'; 'MM Cap'; 'Accelerating'; 'Delaying'; 'Transitory'; 'Final Tariffs'; 'Initial Tariffs'; 'No Bond'; 'Low Bond Adjustment' ;'Low Habit'; 'Two Sectors'; 'NoIO'; 'Common Aggreg'; 'Static'; 'Static FA'};
Static = [Delta_Wel_SOE_Het_ss; Delta_Wel_SOE_sectorAG_ss; Delta_Wel_SOE_sectorMM1_ss; Delta_Wel_SOE_sectorMM3_ss; Delta_Wel_SOE_sectorMM4_ss; Delta_Wel_SOE_HetAcc_ss; Delta_Wel_SOE_HetDel_ss; Delta_Wel_SOE_HetBack_ss; Delta_Wel_SOE_Hom_ss; Delta_Wel_SOE_HetHomIni_ss; Delta_Wel_SOE_HetNoBond_ss; Delta_Wel_SOE_HetNoBAC_ss; Delta_Wel_SOE_HetNoHab_ss; Delta_Wel_SOE_TwoS_ss; Delta_Wel_SOE_HetNoIO_ss; Delta_Wel_SOE_HetIC_ss; Delta_Wel_SOE_HetStatic_ss; Delta_Wel_SOE_HetStaticFA_ss]*100;
Dynamic = [Delta_Wel_SOE_Het;  Delta_Wel_SOE_sectorAG; Delta_Wel_SOE_sectorMM1; Delta_Wel_SOE_sectorMM3; Delta_Wel_SOE_sectorMM4; Delta_Wel_SOE_HetAcc; Delta_Wel_SOE_HetDel; Delta_Wel_SOE_HetBack; Delta_Wel_SOE_Hom; Delta_Wel_SOE_HetHomIni; Delta_Wel_SOE_HetNoBond; Delta_Wel_SOE_HetNoBAC; Delta_Wel_SOE_HetNoHab; Delta_Wel_SOE_TwoS; Delta_Wel_SOE_HetNoIO; Delta_Wel_SOE_HetIC; Delta_Wel_SOE_HetStatic; Delta_Wel_SOE_HetStaticFA]*100;

T = table(Moment,Static, Dynamic)
