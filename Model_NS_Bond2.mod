%% Organize variables
%Aggregates country
%@#define Sectors = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
@#define Sectors = ["1", "2", "3", "4", "5", "6"]
@#define Sectors2 = ["1", "3", "4", "5", "6"]
@#define Country = ["SOE", "ROW"]

%----------------------------------------------------------------
% define variables 
%----------------------------------------------------------------

%% Country
var

R_bnd

@#for co in Country
   P_C_@{co} Q_@{co} C_@{co} W_@{co} NT_@{co} TauT_@{co} EX_T_@{co} PI_T_@{co} 
   WL_@{co} D_C_@{co} NE_@{co} IM_T_@{co} WR_@{co} EX_T_@{co}_real IM_T_@{co}_real
   Lp_T_@{co} GO_T_@{co} I_T_@{co} I_T_@{co}_real B_@{co} lambda_@{co}
@#endfor


%Sectors
@#for co in Country
@#for sec in Sectors

    P_@{co}@{sec} K_@{co}@{sec} P_K_@{co}@{sec} R_K_@{co}@{sec} X_C_@{co}@{sec} I_@{co}@{sec} 
    D_K_@{co}@{sec} X_@{co}@{sec} D_@{co}@{sec} MC_@{co}@{sec} PI_0_@{co}@{sec} a_p_@{co}@{sec} 
    a_0_@{co}@{sec} a_1_@{co}@{sec} Psi_p_@{co}@{sec} Psi_0_@{co}@{sec} Psi_1_@{co}@{sec} 
    n_p_@{co}@{sec} n_0_@{co}@{sec} n_1_@{co}@{sec} EV_0_@{co}@{sec} EV_1_@{co}@{sec}  
    NX_@{co}@{sec} NX_1_@{co}@{sec} NX_0_@{co}@{sec} N_0_@{co}@{sec} Lp_@{co}@{sec} 
    Lf_@{co}@{sec} TauR_@{co}@{sec} EX_@{co}@{sec} PI_@{co}@{sec} 
    EV_inf_@{co}@{sec} mu_@{co}@{sec} IM_@{co}@{sec} 
    GO_@{co}@{sec} P_M_@{co}@{sec} mm_@{co}@{sec} M_@{co}@{sec} D_M_@{co}@{sec}
    MC_CD_@{co}@{sec} CD_@{co}@{sec} EX_@{co}@{sec}_real PI_@{co}@{sec}_real
    GO_dom_@{co}@{sec} EX_0_@{co}@{sec} EX_1_@{co}@{sec} EX_0_@{co}@{sec}_real EX_1_@{co}@{sec}_real
    TO_@{co}@{sec} TB_@{co}@{sec} D_@{co}@{sec}_real IM_@{co}@{sec}_real
@#for sec2 in Sectors 
  X_K_@{co}@{sec}@{sec2}
  X_M_@{co}@{sec}@{sec2}
@#endfor

@#endfor
@#endfor

@#for sec in Sectors
    tau_SOE@{sec}ROW
    tau_ROW@{sec}SOE  
@#endfor

B_SOE_ss B_ROW_ss

;

varexo

@#for sec in Sectors
    eps_SOE@{sec}ROW  
    eps_ROW@{sec}SOE
@#endfor


eps_B_SOE_ss eps_B_ROW_ss
;

%----------------------------------------------------------------
% define parameters
%----------------------------------------------------------------

parameters 

%Aggregate
beta eta theta S

%Country
@#for co in Country  
    L_@{co} n_@{co} fe_@{co} W_@{co}_ss 
@#endfor


%Sectors All
@#for co in Country  
@#for sec in Sectors

    sigma_@{co}@{sec} delta_@{co}@{sec} omega_c_@{co}@{sec} alpha_@{co}@{sec} 
    theta_s_@{co}@{sec} f0_@{co}@{sec} f1_@{co}@{sec} fp_@{co}@{sec} 
    phi_k_@{co}@{sec} kappa_@{co}@{sec} sigma_m_@{co}@{sec} 
    P_@{co}@{sec}_ss
    P_K_@{co}@{sec}_ss
    
%Capital
   
@#for sec2 in Sectors
    omega_k_@{co}@{sec}@{sec2}
    lambda_@{co}@{sec}@{sec2}
@#endfor
       
@#for co2 in Country  
    omega_@{co}@{sec}@{co2}
@#endfor
@#endfor
@#endfor

@#for sec in Sectors
    xi_0_SOE@{sec}ROW
    xi_0_ROW@{sec}SOE
    
    xi_1_SOE@{sec}ROW    
    xi_1_ROW@{sec}SOE 
    
@#endfor

phi_b phi_c phi_bb

;


%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------

load Params_AG;
load Params_CO;
load Params_S;
load Params_I;
load Params_Omega;
load Params_lambda;

load Xss0_CO
load Xss0_S
load Xss0_K
load Xss0_M

load Tau0
load Ice0_0
load Ice1_0
load Ice0_1
load Ice1_1

load Xss1_CO
load Xss1_S
load Xss1_K
load Xss1_M
load Tau1

beta = Params_AG(1,1);
eta  = Params_AG(2,1);
theta = Params_AG(3,1);
S = Params_AG(4,1);

L_SOE = Params_CO(1,1); 
n_SOE = Params_CO(2,1);
fe_SOE = Params_CO(3,1);

W_SOE_ss  = Xss0_CO(4,1);
W_ROW_ss  = Xss0_CO(4,2);


L_ROW = Params_CO(1,2); 
n_ROW = Params_CO(2,2);
fe_ROW = Params_CO(3,2);

@#for sec in Sectors
omega_SOE@{sec}SOE = Params_Omega(1,@{sec},1);
omega_SOE@{sec}ROW = Params_Omega(1,@{sec},2);

omega_ROW@{sec}SOE = Params_Omega(2,@{sec},1);
omega_ROW@{sec}ROW = Params_Omega(2,@{sec},2);

@#endfor

@#for sec in Sectors

    sigma_SOE@{sec} = Params_S(1,1,@{sec}); 
    delta_SOE@{sec} = Params_S(2,1,@{sec});
    omega_c_SOE@{sec} = Params_S(3,1,@{sec});
    alpha_SOE@{sec} = Params_S(4,1,@{sec});
    theta_s_SOE@{sec} = Params_S(5,1,@{sec});
    fp_SOE@{sec} = Params_S(6,1,@{sec});
    f0_SOE@{sec} = Params_S(7,1,@{sec});
    f1_SOE@{sec} = Params_S(8,1,@{sec});
    %phi_k_SOE@{sec} = Params_S(9,1,@{sec});
    kappa_SOE@{sec} = Params_S(10,1,@{sec});
    sigma_m_SOE@{sec} = Params_S(11,1,@{sec});

    sigma_ROW@{sec} = Params_S(1,2,@{sec}); 
    delta_ROW@{sec} = Params_S(2,2,@{sec});
    omega_c_ROW@{sec} = Params_S(3,2,@{sec});
    alpha_ROW@{sec} = Params_S(4,2,@{sec});
    theta_s_ROW@{sec} = Params_S(5,2,@{sec});
    fp_ROW@{sec} = Params_S(6,2,@{sec});
    f0_ROW@{sec} = Params_S(7,2,@{sec});
    f1_ROW@{sec} = Params_S(8,2,@{sec});
%    phi_k_ROW@{sec} = Params_S(9,2,@{sec});
    
    kappa_ROW@{sec} = Params_S(10,2,@{sec});
    sigma_m_ROW@{sec} = Params_S(11,2,@{sec});

@#endfor


@#for sec in Sectors
    phi_k_SOE@{sec} = 1/25;    
    phi_k_ROW@{sec} = 1/25;   
    
@#endfor

@#for sec in Sectors
    P_SOE@{sec}_ss = Xss0_S(1,1,@{sec});
    P_ROW@{sec}_ss = Xss0_S(1,2,@{sec});
@#endfor

@#for sec in Sectors
@#for sec2 in Sectors

    omega_k_SOE@{sec}@{sec2} = Params_I(1,@{sec},@{sec2});
    omega_k_ROW@{sec}@{sec2} = Params_I(2,@{sec},@{sec2});

    lambda_SOE@{sec}@{sec2} = Params_lambda(1,@{sec},@{sec2});
    lambda_ROW@{sec}@{sec2} = Params_lambda(2,@{sec},@{sec2});   
    
@#endfor
@#endfor

@#for sec in Sectors
    P_K_SOE@{sec}_ss = Xss0_S(3,1,@{sec});
    P_K_ROW@{sec}_ss = Xss0_S(3,2,@{sec});
@#endfor

@#for sec in Sectors
xi_0_SOE@{sec}ROW = Ice0_0(1,@{sec},2);
xi_0_ROW@{sec}SOE = Ice0_0(2,@{sec},1);

xi_1_SOE@{sec}ROW = Ice1_0(1,@{sec},2);
xi_1_ROW@{sec}SOE = Ice1_0(2,@{sec},1);
@#endfor

phi_b = 0.1;
phi_bb = 0.25;
phi_c = 0.75;

%phi_b = 0.0001;
%phi_bb = 0.05;
%phi_c = 0.75;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------
model; 

[name='Price of C']
P_C_SOE = 1;

@#for co in Country
WR_@{co} = W_@{co}/P_C_@{co};
@#endfor

@#for co in Country
[name='P_C*C']
D_C_@{co} = P_C_@{co}*C_@{co};
[name='Q']
%Q_@{co} = beta*(C_@{co})/(C_@{co}(+1));
%Q_@{co} = beta*(C_@{co}-phi_c*C_@{co}(-1))/(C_@{co}(+1)-phi_c*C_@{co});
Q_@{co} = beta*(C_@{co}-phi_c*C_@{co}(-1))^(2)/(C_@{co}(+1)-phi_c*C_@{co})^(2);
@#endfor

@#for co in Country
[name='Euler equation Bonds']
%1/C_@{co}/P_C_@{co} = lambda_@{co};
%1/(C_@{co}-phi_c*C_@{co}(-1)) - P_C_@{co}*lambda_@{co} - phi_c*beta/(C_@{co}(+1)-phi_c*C_@{co});
(C_@{co}-phi_c*C_@{co}(-1))^(-2) - P_C_@{co}*lambda_@{co} - phi_c*beta*(C_@{co}(+1)-phi_c*C_@{co})^(-2);

lambda_@{co}*(1+phi_b*(B_@{co}-B_@{co}_ss)+phi_bb*(B_@{co}-B_@{co}(-1))) = beta*lambda_@{co}(+1)*(R_bnd+phi_bb*(B_@{co}(+1)-B_@{co}));


@#for sec in Sectors
[name='Euler equation']
%mu_@{co}@{sec} = beta*(1/C_@{co}(+1)/P_C_@{co}(+1)*R_K_@{co}@{sec}(+1) + mu_@{co}@{sec}(+1)*(1-delta_@{co}@{sec}));
%1/C_@{co}*P_K_@{co}@{sec}/P_C_@{co} = mu_@{co}@{sec}*(1-phi_k_@{co}@{sec}*(I_@{co}@{sec}/I_@{co}@{sec}(-1)-1)/I_@{co}@{sec}(-1)) + beta*mu_@{co}@{sec}(+1)*phi_k_@{co}@{sec}*(I_@{co}@{sec}(+1)/I_@{co}@{sec}-1)*I_@{co}@{sec}(+1)/I_@{co}@{sec}^2;

mu_@{co}@{sec} = beta*(lambda_@{co}(+1)*R_K_@{co}@{sec}(+1) + mu_@{co}@{sec}(+1)*(1-delta_@{co}@{sec}));
%lambda_@{co}*P_K_@{co}@{sec} = mu_@{co}@{sec}*(1-phi_k_@{co}@{sec}*(I_@{co}@{sec}/I_@{co}@{sec}(-1)-1)/I_@{co}@{sec}(-1)) + beta*mu_@{co}@{sec}(+1)*phi_k_@{co}@{sec}*(I_@{co}@{sec}(+1)/I_@{co}@{sec}-1)*I_@{co}@{sec}(+1)/I_@{co}@{sec}^2;
lambda_@{co}*P_K_@{co}@{sec}*(1+phi_k_@{co}@{sec}*(K_@{co}@{sec}-K_@{co}@{sec}(-1))) = beta*lambda_@{co}(+1)*(R_K_@{co}@{sec}(+1)+(1-delta_@{co}@{sec})*P_K_@{co}@{sec}(+1)+phi_k_@{co}@{sec}*(K_@{co}@{sec}(+1)-K_@{co}@{sec})*P_K_@{co}@{sec}(+1));

[name='Demand for sector-interm: C']
X_C_@{co}@{sec} = omega_c_@{co}@{sec}*P_@{co}@{sec}^(-theta)*P_C_@{co}^(theta-1)*D_C_@{co};
[name='Investment']
%I_@{co}@{sec} = K_@{co}@{sec}-(1-delta_@{co}@{sec})*K_@{co}@{sec}(-1) + phi_k_@{co}@{sec}/2*(I_@{co}@{sec}/I_@{co}@{sec}(-1)-1)^2;
I_@{co}@{sec} = K_@{co}@{sec}-(1-delta_@{co}@{sec})*K_@{co}@{sec}(-1) + phi_k_@{co}@{sec}/2*(K_@{co}@{sec}-K_@{co}@{sec}(-1))^2;
[name='P_K*I']
D_K_@{co}@{sec} = I_@{co}@{sec}*P_K_@{co}@{sec};
[name='D_M']
D_M_@{co}@{sec} = M_@{co}@{sec}*P_M_@{co}@{sec};
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Price of K']
P_K_@{co}@{sec}^(1-sigma_@{co}@{sec}) = 
@#for sec2 in Sectors
    +omega_k_@{co}@{sec}@{sec2}*P_@{co}@{sec2}^(1-sigma_@{co}@{sec})
@#endfor
;
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Price of M']
P_M_@{co}@{sec}^(1-sigma_m_@{co}@{sec}) = 
@#for sec2 in Sectors
    +lambda_@{co}@{sec}@{sec2}*P_@{co}@{sec2}^(1-sigma_m_@{co}@{sec})
@#endfor
;
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
@#for sec2 in Sectors
[name='Intermediates K']
    X_K_@{co}@{sec}@{sec2} = omega_k_@{co}@{sec}@{sec2}*P_@{co}@{sec2}^(-sigma_@{co}@{sec})*P_K_@{co}@{sec}^(sigma_@{co}@{sec}-1)*D_K_@{co}@{sec};   
[name='Intermediates M']
    X_M_@{co}@{sec}@{sec2} = lambda_@{co}@{sec}@{sec2}*P_@{co}@{sec2}^(-sigma_m_@{co}@{sec})*P_M_@{co}@{sec}^(sigma_m_@{co}@{sec}-1)*D_M_@{co}@{sec};   
@#endfor
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Total demand for intermediates']
    X_@{co}@{sec} = X_C_@{co}@{sec} 
@#for sec2 in Sectors
    + X_K_@{co}@{sec2}@{sec} + X_M_@{co}@{sec2}@{sec}
@#endfor
;
[name='Nominal demand for intermediates']
D_@{co}@{sec} = P_@{co}@{sec}*X_@{co}@{sec};
D_@{co}@{sec}_real = P_@{co}@{sec}_ss*X_@{co}@{sec};
[name='Marginal Cost CD']
MC_CD_@{co}@{sec} = (R_K_@{co}@{sec}/alpha_@{co}@{sec})^alpha_@{co}@{sec}*(W_@{co}/kappa_@{co}@{sec})^kappa_@{co}@{sec}*(P_M_@{co}@{sec}/(1-alpha_@{co}@{sec}-kappa_@{co}@{sec}))^(1-alpha_@{co}@{sec}-kappa_@{co}@{sec});
[name='Marginal Cost']
MC_@{co}@{sec} = MC_CD_@{co}@{sec};
[name='PI 0']
PI_0_@{co}@{sec} = 1/theta_s_@{co}@{sec}*(theta_s_@{co}@{sec}*MC_@{co}@{sec}/(theta_s_@{co}@{sec}-1))^(1-theta_s_@{co}@{sec});
[name='Marginal producer']
a_p_@{co}@{sec} = W_@{co}*(fp_@{co}@{sec})/(PI_0_@{co}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co});

@#endfor
@#endfor



@#for sec in Sectors
@#for co in Country
[name='Marginal new exporter']
W_@{co}*f0_@{co}@{sec} = n_@{co}*Q_@{co}*(EV_1_@{co}@{sec}(+1)-EV_inf_@{co}@{sec}(+1))*(NE_@{co}/(NE_@{co}(-1)+NE_@{co}(-2))*2)^0
@#for co1 in Country
    @#if co1!=co
    + PI_0_@{co}@{sec}*a_0_@{co}@{sec}*(xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})/(tau_@{co1}@{sec}@{co}^theta_s_@{co1}@{sec}*P_@{co1}@{sec}^(1-theta_s_@{co1}@{sec})))
    @#endif
@#endfor
;
@#endfor

@#for co in Country
[name='Marginal continuation exporter']
W_@{co}*f1_@{co}@{sec} = n_@{co}*Q_@{co}*(EV_1_@{co}@{sec}(+1)-EV_inf_@{co}@{sec}(+1))
@#for co1 in Country
    @#if co1!=co
    + PI_0_@{co}@{sec}*a_1_@{co}@{sec}*(xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})/(tau_@{co1}@{sec}@{co}^theta_s_@{co1}@{sec}*P_@{co1}@{sec}^(1-theta_s_@{co1}@{sec})))
    @#endif
@#endfor
;
@#endfor

@#for co in Country
[name='EV_0']
EV_0_@{co}@{sec} = PI_0_@{co}@{sec}*(P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co}*Psi_p_@{co}@{sec}) - W_@{co}*(fp_@{co}@{sec})*n_p_@{co}@{sec} - W_@{co}*f0_@{co}@{sec}*n_0_@{co}@{sec} + n_@{co}*Q_@{co}*((1-n_1_@{co}@{sec}(+1))*EV_inf_@{co}@{sec}(+1) + n_1_@{co}@{sec}(+1)*EV_1_@{co}@{sec}(+1))
@#for co1 in Country
    @#if co1!=co
    + PI_0_@{co}@{sec}*(xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})*Psi_0_@{co}@{sec}) 
    @#endif
@#endfor
;
@#endfor


@#for co in Country
[name='EV_1']
EV_1_@{co}@{sec} = PI_0_@{co}@{sec}*(P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co}*Psi_p_@{co}@{sec}) - W_@{co}*(fp_@{co}@{sec})*n_p_@{co}@{sec} - W_@{co}*f1_@{co}@{sec}*n_1_@{co}@{sec} + n_@{co}*Q_@{co}*((1-n_1_@{co}@{sec}(+1))*EV_inf_@{co}@{sec}(+1) + n_1_@{co}@{sec}(+1)*EV_1_@{co}@{sec}(+1))
@#for co1 in Country
    @#if co1!=co
    + PI_0_@{co}@{sec}*(xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})*Psi_1_@{co}@{sec}) 
    @#endif
@#endfor
;
@#endfor
@#endfor

                                                                                                 
@#for sec in Sectors
@#for co in Country                                                                                                                        
[name='EV_inf']
EV_inf_@{co}@{sec} = PI_0_@{co}@{sec}*(P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co}*Psi_p_@{co}@{sec})
                - W_@{co}*(fp_@{co}@{sec})*n_p_@{co}@{sec} + n_@{co}*Q_@{co}*((1-n_0_@{co}@{sec}(+1))*EV_inf_@{co}@{sec}(+1) + n_0_@{co}@{sec}(+1)*EV_0_@{co}@{sec}(+1));               
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors

[name='Psi_p']
Psi_p_@{co}@{sec} = eta/(eta-1)*a_p_@{co}@{sec}^(1-eta);
[name='Psi_0']
Psi_0_@{co}@{sec} = eta/(eta-1)*a_0_@{co}@{sec}^(1-eta);
[name='Psi_1']
Psi_1_@{co}@{sec} = eta/(eta-1)*a_1_@{co}@{sec}^(1-eta);

[name='n_p']
n_p_@{co}@{sec} = a_p_@{co}@{sec}^(-eta);
[name='n_0']
n_0_@{co}@{sec} = a_0_@{co}@{sec}^(-eta);
[name='n_1']
n_1_@{co}@{sec} = a_1_@{co}@{sec}^(-eta);
@#endfor
@#endfor

[name='Free entry']
@#for co in Country
W_@{co}*fe_@{co} =
@#for sec in Sectors
 %   + Q_@{co}*EV_inf_@{co}@{sec}(+1)/S/(NT_@{co}(+1)/(NT_@{co}+NT_@{co}(-1))*2)^0
   + Q_@{co}*EV_inf_@{co}@{sec}(+1)/S/(NE_@{co}/(NE_@{co}(-1)+NE_@{co}(-2))*2)
@#endfor
;
@#endfor

@#for co in Country
[name='Masses of Firms']
    NT_@{co} = n_@{co}*NT_@{co}(-1) + NE_@{co}(-1);
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Masses of Exporters']
    NX_@{co}@{sec} = NX_0_@{co}@{sec} + NX_1_@{co}@{sec};
[name='Masses of Exporters (Cont)']
    NX_1_@{co}@{sec} = n_@{co}*n_1_@{co}@{sec}*NX_@{co}@{sec}(-1)*(1+0.20*(NT_@{co}(+1)-NT_@{co}));
[name='Masses of Exporters (New)']
    NX_0_@{co}@{sec}/(1+0.20*(NT_@{co}(+1)-NT_@{co})) = n_0_@{co}@{sec}*(NE_@{co}(-1)/S+n_@{co}*N_0_@{co}@{sec}(-1));
[name='Masses of Non-Exporters']
    NT_@{co}/S = N_0_@{co}@{sec} + NX_@{co}@{sec};
@#endfor
@#endfor

%%%%%%%%%%%%%%%%%%%%
@#for co in Country
@#for sec in Sectors
[name='Production Labor']
Lp_@{co}@{sec} = (theta_s_@{co}@{sec}-1)*PI_0_@{co}@{sec}*(kappa_@{co}@{sec})/W_@{co}*(NT_@{co}/S*Psi_p_@{co}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co})
@#for co1 in Country
    @#if co1!=co
          + (theta_s_@{co}@{sec}-1)*PI_0_@{co}@{sec}*(kappa_@{co}@{sec})/W_@{co}*(tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})*(xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+ 
           xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec}))
    @#endif
@#endfor
;
@#endfor
@#endfor
%%%%%%%%%%%%%%%%%%%%

@#for co in Country
@#for sec in Sectors
[name='Labor for fixed costs']
Lf_@{co}@{sec} = NT_@{co}/S*n_p_@{co}@{sec}*(fp_@{co}@{sec}) + f0_@{co}@{sec}*NX_0_@{co}@{sec} + f1_@{co}@{sec}*NX_1_@{co}@{sec} + NE_@{co}*fe_@{co}/S;
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Tariff Rev']
TauR_@{co}@{sec} = 
@#for co1 in Country
    @#if co1!=co
    +(omega_@{co}@{sec}@{co1})*(tau_@{co}@{sec}@{co1}-1)*theta_s_@{co}@{sec}*PI_0_@{co1}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*tau_@{co}@{sec}@{co1}^(-theta_s_@{co}@{sec})*
            (xi_1_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_1_@{co1}@{sec}/n_1_@{co1}@{sec}*Psi_1_@{co1}@{sec}+ 
             xi_0_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_0_@{co1}@{sec}/n_0_@{co1}@{sec}*Psi_0_@{co1}@{sec})          
    @#endif
@#endfor
;
@#endfor
@#endfor

@#for co in Country
@#for sec in Sectors
[name = 'Imports']             
IM_@{co}@{sec} =  
@#for co1 in Country
    @#if co1!=co            
    + omega_@{co}@{sec}@{co1}*theta_s_@{co}@{sec}*PI_0_@{co1}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*tau_@{co}@{sec}@{co1}^(-theta_s_@{co}@{sec})*
            (xi_1_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_1_@{co1}@{sec}/n_1_@{co1}@{sec}*Psi_1_@{co1}@{sec}+
             xi_0_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_0_@{co1}@{sec}/n_0_@{co1}@{sec}*Psi_0_@{co1}@{sec})
    @#endif
@#endfor
;
@#endfor
@#endfor
             

@#for co in Country
@#for sec in Sectors
[name = 'Imports']             
IM_@{co}@{sec}_real =  
@#for co1 in Country
    @#if co1!=co            
    + omega_@{co}@{sec}@{co1}*theta_s_@{co}@{sec}*PI_0_@{co1}@{sec}*P_@{co}@{sec}_ss^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}_real*tau_@{co}@{sec}@{co1}^(-theta_s_@{co}@{sec})*
            (xi_1_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_1_@{co1}@{sec}/n_1_@{co1}@{sec}*Psi_1_@{co1}@{sec}+
             xi_0_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_0_@{co1}@{sec}/n_0_@{co1}@{sec}*Psi_0_@{co1}@{sec})
    @#endif
@#endfor
;
@#endfor
@#endfor


@#for co in Country
@#for sec in Sectors
[name = 'Gross Output']             
GO_@{co}@{sec} =  theta_s_@{co}@{sec}*PI_0_@{co}@{sec}*(NT_@{co}/S*Psi_p_@{co}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co})
@#for co1 in Country
    @#if co1!=co            
        + theta_s_@{co}@{sec}*PI_0_@{co}@{sec}*(tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})*
          (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+
           xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec}))
    @#endif
@#endfor
;
@#endfor
@#endfor


@#for sec in Sectors            
[name='Exports']
@#for co in Country
EX_@{co}@{sec} =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+
             xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor


[name='Exports Old']
@#for co in Country
EX_1_@{co}@{sec} =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor

[name='Exports New']
@#for co in Country
EX_0_@{co}@{sec} =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor


[name='Exports real']
@#for co in Country
EX_@{co}@{sec}_real =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}_ss^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}_real*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+
             xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor

[name='Exports real (old)']
@#for co in Country
EX_1_@{co}@{sec}_real =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}_ss^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}_real*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor

[name='Exports real (new)']
@#for co in Country
EX_0_@{co}@{sec}_real =
@#for co1 in Country
    @#if co1!=co
           +(omega_@{co1}@{sec}@{co})*theta_s_@{co1}@{sec}*PI_0_@{co}@{sec}*P_@{co1}@{sec}_ss^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}_real*tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*
            (xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co1}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec})
    @#endif
@#endfor
;
@#endfor          

@#for co in Country
GO_dom_@{co}@{sec} = theta_s_@{co}@{sec}*PI_0_@{co}@{sec}*(NT_@{co}/S*Psi_p_@{co}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co});
@#endfor           
             
[name='Profits']
@#for co in Country
PI_@{co}@{sec} = PI_0_@{co}@{sec}*(NT_@{co}/S*Psi_p_@{co}@{sec}*P_@{co}@{sec}^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}*omega_@{co}@{sec}@{co})- W_@{co}*Lf_@{co}@{sec} 
@#for co1 in Country
    @#if co1!=co
           + PI_0_@{co}@{sec}*(tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}*(omega_@{co1}@{sec}@{co})*
          (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+
           xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec}))
    @#endif
@#endfor
;
@#endfor


[name='Profits real']
@#for co in Country
PI_@{co}@{sec}_real = PI_0_@{co}@{sec}*(NT_@{co}/S*Psi_p_@{co}@{sec}*P_@{co}@{sec}_ss^(theta_s_@{co}@{sec}-1)*D_@{co}@{sec}_real*omega_@{co}@{sec}@{co})- W_@{co}_ss*Lf_@{co}@{sec} 
@#for co1 in Country
    @#if co1!=co
           + PI_0_@{co}@{sec}*(tau_@{co1}@{sec}@{co}^(-theta_s_@{co1}@{sec})*P_@{co1}@{sec}_ss^(theta_s_@{co1}@{sec}-1)*D_@{co1}@{sec}_real*(omega_@{co1}@{sec}@{co})*
          (xi_1_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_1_@{co}@{sec}/n_1_@{co}@{sec}*Psi_1_@{co}@{sec}+
           xi_0_@{co}@{sec}@{co1}^(1-theta_s_@{co}@{sec})*NX_0_@{co}@{sec}/n_0_@{co}@{sec}*Psi_0_@{co}@{sec}))
    @#endif
@#endfor
;
@#endfor            
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Sector Price'] 
P_@{co}@{sec}^(1-theta_s_@{co}@{sec}) = omega_@{co}@{sec}@{co}*NT_@{co}/S*(theta_s_@{co}@{sec}*MC_@{co}@{sec}/(theta_s_@{co}@{sec}-1))^(1-theta_s_@{co}@{sec})*Psi_p_@{co}@{sec} 
@#for co1 in Country
    @#if co1!=co      
      +  omega_@{co}@{sec}@{co1}*(tau_@{co}@{sec}@{co1}*theta_s_@{co}@{sec}/(theta_s_@{co}@{sec}-1)*MC_@{co1}@{sec})^(1-theta_s_@{co}@{sec})*
       (xi_1_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_1_@{co1}@{sec}/n_1_@{co1}@{sec}*Psi_1_@{co1}@{sec}+ 
        xi_0_@{co1}@{sec}@{co}^(1-theta_s_@{co}@{sec})*NX_0_@{co1}@{sec}/n_0_@{co1}@{sec}*Psi_0_@{co1}@{sec})
    @#endif
@#endfor
;
@#endfor            
@#endfor


@#for co in Country
[name='Tariffs Total']
    TauT_@{co} = 0
@#for sec in Sectors
    + TauR_@{co}@{sec}
@#endfor
;
@#endfor

@#for co in Country
[name='Exports Total']
    EX_T_@{co} = 0
@#for sec in Sectors
    + EX_@{co}@{sec}
@#endfor
;
@#endfor

@#for co in Country
[name='Lp Total']
Lp_T_@{co} = 0 
@#for sec in Sectors
    + Lp_@{co}@{sec}
@#endfor
;
@#endfor


@#for co in Country
[name='GO Total']
GO_T_@{co} = 0 
@#for sec in Sectors
    + GO_@{co}@{sec}
@#endfor
;
@#endfor


@#for co in Country
[name='Exports Total']
    EX_T_@{co}_real = 0
@#for sec in Sectors
    + EX_@{co}@{sec}_real
@#endfor
;
@#endfor

@#for co in Country
[name='Imports Total']
    IM_T_@{co} = 0
@#for sec in Sectors
    + IM_@{co}@{sec}
@#endfor
;
@#endfor

@#for co in Country
[name='Imports Total real']
    IM_T_@{co}_real = 0
@#for sec in Sectors
    + IM_@{co}@{sec}_real
@#endfor
;
@#endfor


@#for co in Country
[name='Profits Total']
    PI_T_@{co} = 0
@#for sec in Sectors
    + PI_@{co}@{sec}
@#endfor
;
@#endfor


@#for co in Country
[name='I_Total']
I_T_@{co} = 0
@#for sec in Sectors
    + P_K_@{co}@{sec}*I_@{co}@{sec}
@#endfor
;
@#endfor

@#for co in Country
[name='I_Total real']
I_T_@{co}_real = 0
@#for sec in Sectors
    + P_K_@{co}@{sec}_ss*I_@{co}@{sec}
@#endfor
;
@#endfor

@#for co in Country
[name='Walras']
    WL_@{co} = P_C_@{co}*C_@{co} + B_@{co} + phi_b/2*(B_@{co}-B_@{co}_ss)^2 + phi_bb/2*(B_@{co}-B_@{co}(-1))^2 - R_bnd*B_@{co}(+1) - TauT_@{co} - PI_T_@{co} - W_@{co}*L_@{co}  
@#for sec in Sectors
     + P_K_@{co}@{sec}*I_@{co}@{sec} - R_K_@{co}@{sec}*K_@{co}@{sec}(-1)       
@#endfor
;
@#endfor


@#for co in Country
[name='Total Labor']        
L_@{co} = 0
@#for sec in Sectors
  + Lp_@{co}@{sec} + Lf_@{co}@{sec}  
@#endfor
;
@#endfor

@#for co in Country
@#for sec in Sectors
[name='Capital']        
K_@{co}@{sec}(-1) = alpha_@{co}@{sec}/(kappa_@{co}@{sec})*W_@{co}/R_K_@{co}@{sec}*Lp_@{co}@{sec};
[name='mm']        
mm_@{co}@{sec} = (1-kappa_@{co}@{sec}-alpha_@{co}@{sec})/(kappa_@{co}@{sec})*W_@{co}/P_M_@{co}@{sec}*Lp_@{co}@{sec};
[name= 'M']
M_@{co}@{sec} = mm_@{co}@{sec};
[name='CD']
CD_@{co}@{sec} = K_@{co}@{sec}(-1)^alpha_@{co}@{sec}*Lp_@{co}@{sec}^kappa_@{co}@{sec}*mm_@{co}@{sec}^(1-alpha_@{co}@{sec}-kappa_@{co}@{sec});
[name='TO']
TO_@{co}@{sec} = (EX_@{co}@{sec}+IM_@{co}@{sec})/GO_@{co}@{sec};
[name='TB']
TB_@{co}@{sec} = (EX_@{co}@{sec}-IM_@{co}@{sec})/GO_@{co}@{sec};
@#endfor
@#endfor

[name= 'P_C']
@#for co in Country
P_C_@{co}^(1-theta) = 
@#for sec in Sectors
   + omega_c_@{co}@{sec}*P_@{co}@{sec}^(1-theta)
@#endfor
;
@#endfor

[name= 'Trade Balance']
%EX_T_SOE = IM_T_SOE;
%WL_SOE = 0;
%WL_ROW = 0;
EX_T_SOE - IM_T_SOE = R_bnd*B_SOE(+1) - (B_SOE + phi_b/2*(B_SOE-B_SOE_ss)^2 + phi_bb/2*(B_SOE-B_SOE(-1))^2);
B_SOE + B_ROW = 0;

[name= 'Tariffs']
@#for sec in Sectors
    tau_SOE@{sec}ROW = eps_SOE@{sec}ROW;
    tau_ROW@{sec}SOE = eps_ROW@{sec}SOE;
@#endfor

B_SOE_ss = eps_B_SOE_ss;
B_ROW_ss = eps_B_ROW_ss;

end;

%----------------------------------------------------------------
%  Initial values
%---------------------------------------------------------------


initval;

R_bnd = 1/beta;

P_C_SOE = Xss0_CO(1,1); 
Q_SOE  = Xss0_CO(2,1);
C_SOE  = Xss0_CO(3,1);
W_SOE  = Xss0_CO(4,1);
NT_SOE = Xss0_CO(5,1); 
TauT_SOE = Xss0_CO(6,1);
EX_T_SOE = Xss0_CO(7,1);
PI_T_SOE = Xss0_CO(8,1);
WL_SOE = Xss0_CO(9,1);
D_C_SOE = Xss0_CO(10,1);
NE_SOE = Xss0_CO(11,1); 
IM_T_SOE = Xss0_CO(12,1);
Lp_T_SOE = Xss0_CO(13,1);

lambda_SOE = 1/C_SOE/P_C_SOE;

P_C_ROW = Xss0_CO(1,2); 
Q_ROW  = Xss0_CO(2,2);
C_ROW  = Xss0_CO(3,2);
W_ROW  = Xss0_CO(4,2);
NT_ROW = Xss0_CO(5,2); 
TauT_ROW = Xss0_CO(6,2);
EX_T_ROW = Xss0_CO(7,2);
PI_T_ROW = Xss0_CO(8,2);
WL_ROW = Xss0_CO(9,2);
D_C_ROW = Xss0_CO(10,2);
NE_ROW = Xss0_CO(11,2);
IM_T_ROW = Xss0_CO(12,2); 
Lp_T_ROW = Xss0_CO(13,2);

lambda_ROW = 1/C_ROW/P_C_ROW;

@#for sec in Sectors

    P_SOE@{sec} = Xss0_S(1,1,@{sec});
    K_SOE@{sec} = Xss0_S(2,1,@{sec});
    P_K_SOE@{sec} = Xss0_S(3,1,@{sec});
    R_K_SOE@{sec} = Xss0_S(4,1,@{sec});
    X_C_SOE@{sec} = Xss0_S(5,1,@{sec});
    I_SOE@{sec} = Xss0_S(6,1,@{sec});
    D_K_SOE@{sec} = Xss0_S(7,1,@{sec});
    X_SOE@{sec} = Xss0_S(8,1,@{sec});
    D_SOE@{sec} = Xss0_S(9,1,@{sec});
    MC_SOE@{sec} = Xss0_S(10,1,@{sec});
    PI_0_SOE@{sec} = Xss0_S(11,1,@{sec});
    a_p_SOE@{sec} = Xss0_S(12,1,@{sec});
    a_0_SOE@{sec} = Xss0_S(13,1,@{sec});
    a_1_SOE@{sec} = Xss0_S(14,1,@{sec});
    Psi_p_SOE@{sec} = Xss0_S(15,1,@{sec});
    Psi_0_SOE@{sec} = Xss0_S(16,1,@{sec});
    Psi_1_SOE@{sec} = Xss0_S(17,1,@{sec});
    n_p_SOE@{sec} = Xss0_S(18,1,@{sec});
    n_0_SOE@{sec} = Xss0_S(19,1,@{sec});
    n_1_SOE@{sec} = Xss0_S(20,1,@{sec});
    EV_0_SOE@{sec} = Xss0_S(21,1,@{sec});
    EV_1_SOE@{sec} = Xss0_S(22,1,@{sec});
    NX_SOE@{sec} = Xss0_S(23,1,@{sec});
    NX_1_SOE@{sec} = Xss0_S(24,1,@{sec});
    NX_0_SOE@{sec} = Xss0_S(25,1,@{sec});
    N_0_SOE@{sec} = Xss0_S(26,1,@{sec});
    Lp_SOE@{sec} = Xss0_S(27,1,@{sec});
    Lf_SOE@{sec} = Xss0_S(28,1,@{sec});
    TauR_SOE@{sec} = Xss0_S(29,1,@{sec});
    EX_SOE@{sec} = Xss0_S(30,1,@{sec});
    PI_SOE@{sec} = Xss0_S(31,1,@{sec});
    EV_inf_SOE@{sec} = Xss0_S(32,1,@{sec});
    mu_SOE@{sec} = Xss0_S(33,1,@{sec});
    IM_SOE@{sec} = Xss0_S(37,1,@{sec});
    P_M_SOE@{sec} =  Xss0_S(38,1,@{sec});
    mm_SOE@{sec} =  Xss0_S(39,1,@{sec});
    M_SOE@{sec} =  Xss0_S(40,1,@{sec});
    D_M_SOE@{sec} =  Xss0_S(41,1,@{sec});
    MC_CD_SOE@{sec} =  Xss0_S(42,1,@{sec});
    CD_SOE@{sec} =  Xss0_S(43,1,@{sec});
    GO_SOE@{sec} =  Xss0_S(44,1,@{sec});

    P_ROW@{sec} = Xss0_S(1,2,@{sec});
    K_ROW@{sec} = Xss0_S(2,2,@{sec});
    P_K_ROW@{sec} = Xss0_S(3,2,@{sec});
    R_K_ROW@{sec} = Xss0_S(4,2,@{sec});
    X_C_ROW@{sec} = Xss0_S(5,2,@{sec});
    I_ROW@{sec} = Xss0_S(6,2,@{sec});
    D_K_ROW@{sec} = Xss0_S(7,2,@{sec});
    X_ROW@{sec} = Xss0_S(8,2,@{sec});
    D_ROW@{sec} = Xss0_S(9,2,@{sec});
    MC_ROW@{sec} = Xss0_S(10,2,@{sec});
    PI_0_ROW@{sec} = Xss0_S(11,2,@{sec});
    a_p_ROW@{sec} = Xss0_S(12,2,@{sec});
    a_0_ROW@{sec} = Xss0_S(13,2,@{sec});
    a_1_ROW@{sec} = Xss0_S(14,2,@{sec});
    Psi_p_ROW@{sec} = Xss0_S(15,2,@{sec});
    Psi_0_ROW@{sec} = Xss0_S(16,2,@{sec});
    Psi_1_ROW@{sec} = Xss0_S(17,2,@{sec});
    n_p_ROW@{sec} = Xss0_S(18,2,@{sec});
    n_0_ROW@{sec} = Xss0_S(19,2,@{sec});
    n_1_ROW@{sec} = Xss0_S(20,2,@{sec});
    EV_0_ROW@{sec} = Xss0_S(21,2,@{sec});
    EV_1_ROW@{sec} = Xss0_S(22,2,@{sec});
    NX_ROW@{sec} = Xss0_S(23,2,@{sec});
    NX_1_ROW@{sec} = Xss0_S(24,2,@{sec});
    NX_0_ROW@{sec} = Xss0_S(25,2,@{sec});
    N_0_ROW@{sec} = Xss0_S(26,2,@{sec});
    Lp_ROW@{sec} = Xss0_S(27,2,@{sec});
    Lf_ROW@{sec} = Xss0_S(28,2,@{sec});
    TauR_ROW@{sec} = Xss0_S(29,2,@{sec});
    EX_ROW@{sec} = Xss0_S(30,2,@{sec});
    PI_ROW@{sec} = Xss0_S(31,2,@{sec});    
    EV_inf_ROW@{sec} = Xss0_S(32,2,@{sec});
    mu_ROW@{sec} = Xss0_S(33,2,@{sec});
    IM_ROW@{sec} = Xss0_S(37,2,@{sec});
    P_M_ROW@{sec} =  Xss0_S(38,2,@{sec});
    mm_ROW@{sec} =  Xss0_S(39,2,@{sec});
    M_ROW@{sec} =  Xss0_S(40,2,@{sec});
    D_M_ROW@{sec} =  Xss0_S(41,2,@{sec});
    MC_CD_ROW@{sec} =  Xss0_S(42,2,@{sec});
    CD_ROW@{sec} =  Xss0_S(43,2,@{sec});
    GO_ROW@{sec} =  Xss0_S(44,2,@{sec});   


@#for sec2 in Sectors
  X_K_SOE@{sec}@{sec2} = Xss0_K(1,@{sec},@{sec2});
  X_K_ROW@{sec}@{sec2} = Xss0_K(2,@{sec},@{sec2});
    
  X_M_SOE@{sec}@{sec2} = Xss0_M(1,@{sec},@{sec2});
  X_M_ROW@{sec}@{sec2} = Xss0_M(2,@{sec},@{sec2});
  
@#endfor


tau_SOE@{sec}ROW = Tau0(1,@{sec},2);
tau_ROW@{sec}SOE = Tau0(2,@{sec},1);

eps_SOE@{sec}ROW = tau_SOE@{sec}ROW; 
eps_ROW@{sec}SOE = tau_ROW@{sec}SOE;

@#endfor

end;
resid;
steady; 
check;
%model_diagnostics;

%%%Final Values

endval;

P_C_SOE = Xss1_CO(1,1); 
Q_SOE  = Xss1_CO(2,1);
C_SOE  = Xss1_CO(3,1);
W_SOE  = Xss1_CO(4,1);
NT_SOE = Xss1_CO(5,1); 
TauT_SOE = Xss1_CO(6,1);
EX_T_SOE = Xss1_CO(7,1);
PI_T_SOE = Xss1_CO(8,1);
WL_SOE = Xss1_CO(9,1);
D_C_SOE = Xss1_CO(10,1);
NE_SOE = Xss1_CO(11,1); 
IM_T_SOE = Xss1_CO(12,1);
Lp_T_SOE = Xss1_CO(13,1);

lambda_SOE = 1/C_SOE/P_C_SOE;

P_C_ROW = Xss1_CO(1,2); 
Q_ROW  = Xss1_CO(2,2);
C_ROW  = Xss1_CO(3,2);
W_ROW  = Xss1_CO(4,2);
NT_ROW = Xss1_CO(5,2); 
TauT_ROW = Xss1_CO(6,2);
EX_T_ROW = Xss1_CO(7,2);
PI_T_ROW = Xss1_CO(8,2);
WL_ROW = Xss1_CO(9,2);
D_C_ROW = Xss1_CO(10,2);
NE_ROW = Xss1_CO(11,2);
IM_T_ROW = Xss1_CO(12,2); 
Lp_T_ROW = Xss1_CO(13,2);

lambda_ROW = 1/C_ROW/P_C_ROW;

@#for sec in Sectors

    P_SOE@{sec} = Xss1_S(1,1,@{sec});
    K_SOE@{sec} = Xss1_S(2,1,@{sec});
    P_K_SOE@{sec} = Xss1_S(3,1,@{sec});
    R_K_SOE@{sec} = Xss1_S(4,1,@{sec});
    X_C_SOE@{sec} = Xss1_S(5,1,@{sec});
    I_SOE@{sec} = Xss1_S(6,1,@{sec});
    D_K_SOE@{sec} = Xss1_S(7,1,@{sec});
    X_SOE@{sec} = Xss1_S(8,1,@{sec});
    D_SOE@{sec} = Xss1_S(9,1,@{sec});
    MC_SOE@{sec} = Xss1_S(10,1,@{sec});
    PI_0_SOE@{sec} = Xss1_S(11,1,@{sec});
    a_p_SOE@{sec} = Xss1_S(12,1,@{sec});
    a_0_SOE@{sec} = Xss1_S(13,1,@{sec});
    a_1_SOE@{sec} = Xss1_S(14,1,@{sec});
    Psi_p_SOE@{sec} = Xss1_S(15,1,@{sec});
    Psi_0_SOE@{sec} = Xss1_S(16,1,@{sec});
    Psi_1_SOE@{sec} = Xss1_S(17,1,@{sec});
    n_p_SOE@{sec} = Xss1_S(18,1,@{sec});
    n_0_SOE@{sec} = Xss1_S(19,1,@{sec});
    n_1_SOE@{sec} = Xss1_S(20,1,@{sec});
    EV_0_SOE@{sec} = Xss1_S(21,1,@{sec});
    EV_1_SOE@{sec} = Xss1_S(22,1,@{sec});
    NX_SOE@{sec} = Xss1_S(23,1,@{sec});
    NX_1_SOE@{sec} = Xss1_S(24,1,@{sec});
    NX_0_SOE@{sec} = Xss1_S(25,1,@{sec});
    N_0_SOE@{sec} = Xss1_S(26,1,@{sec});
    Lp_SOE@{sec} = Xss1_S(27,1,@{sec});
    Lf_SOE@{sec} = Xss1_S(28,1,@{sec});
    TauR_SOE@{sec} = Xss1_S(29,1,@{sec});
    EX_SOE@{sec} = Xss1_S(30,1,@{sec});
    PI_SOE@{sec} = Xss1_S(31,1,@{sec});
    EV_inf_SOE@{sec} = Xss1_S(32,1,@{sec});
    mu_SOE@{sec} = Xss1_S(33,1,@{sec});
    IM_SOE@{sec} = Xss1_S(37,1,@{sec});
    P_M_SOE@{sec} =  Xss1_S(38,1,@{sec});
    mm_SOE@{sec} =  Xss1_S(39,1,@{sec});
    M_SOE@{sec} =  Xss1_S(40,1,@{sec});
    D_M_SOE@{sec} =  Xss1_S(41,1,@{sec});
    MC_CD_SOE@{sec} =  Xss1_S(42,1,@{sec});
    CD_SOE@{sec} =  Xss1_S(43,1,@{sec});
    GO_SOE@{sec} =  Xss1_S(44,1,@{sec});

    P_ROW@{sec} = Xss1_S(1,2,@{sec});
    K_ROW@{sec} = Xss1_S(2,2,@{sec});
    P_K_ROW@{sec} = Xss1_S(3,2,@{sec});
    R_K_ROW@{sec} = Xss1_S(4,2,@{sec});
    X_C_ROW@{sec} = Xss1_S(5,2,@{sec});
    I_ROW@{sec} = Xss1_S(6,2,@{sec});
    D_K_ROW@{sec} = Xss1_S(7,2,@{sec});
    X_ROW@{sec} = Xss1_S(8,2,@{sec});
    D_ROW@{sec} = Xss1_S(9,2,@{sec});
    MC_ROW@{sec} = Xss1_S(10,2,@{sec});
    PI_0_ROW@{sec} = Xss1_S(11,2,@{sec});
    a_p_ROW@{sec} = Xss1_S(12,2,@{sec});
    a_0_ROW@{sec} = Xss1_S(13,2,@{sec});
    a_1_ROW@{sec} = Xss1_S(14,2,@{sec});
    Psi_p_ROW@{sec} = Xss1_S(15,2,@{sec});
    Psi_0_ROW@{sec} = Xss1_S(16,2,@{sec});
    Psi_1_ROW@{sec} = Xss1_S(17,2,@{sec});
    n_p_ROW@{sec} = Xss1_S(18,2,@{sec});
    n_0_ROW@{sec} = Xss1_S(19,2,@{sec});
    n_1_ROW@{sec} = Xss1_S(20,2,@{sec});
    EV_0_ROW@{sec} = Xss1_S(21,2,@{sec});
    EV_1_ROW@{sec} = Xss1_S(22,2,@{sec});
    NX_ROW@{sec} = Xss1_S(23,2,@{sec});
    NX_1_ROW@{sec} = Xss1_S(24,2,@{sec});
    NX_0_ROW@{sec} = Xss1_S(25,2,@{sec});
    N_0_ROW@{sec} = Xss1_S(26,2,@{sec});
    Lp_ROW@{sec} = Xss1_S(27,2,@{sec});
    Lf_ROW@{sec} = Xss1_S(28,2,@{sec});
    TauR_ROW@{sec} = Xss1_S(29,2,@{sec});
    EX_ROW@{sec} = Xss1_S(30,2,@{sec});
    PI_ROW@{sec} = Xss1_S(31,2,@{sec});    
    EV_inf_ROW@{sec} = Xss1_S(32,2,@{sec});
    mu_ROW@{sec} = Xss1_S(33,2,@{sec});
    IM_ROW@{sec} = Xss1_S(37,2,@{sec});
    P_M_ROW@{sec} =  Xss1_S(38,2,@{sec});
    mm_ROW@{sec} =  Xss1_S(39,2,@{sec});
    M_ROW@{sec} =  Xss1_S(40,2,@{sec});
    D_M_ROW@{sec} =  Xss1_S(41,2,@{sec});
    MC_CD_ROW@{sec} =  Xss1_S(42,2,@{sec});
    CD_ROW@{sec} =  Xss1_S(43,2,@{sec});
    GO_ROW@{sec} =  Xss1_S(44,2,@{sec});   


@#for sec2 in Sectors
  X_K_SOE@{sec}@{sec2} = Xss1_K(1,@{sec},@{sec2});
  X_K_ROW@{sec}@{sec2} = Xss1_K(2,@{sec},@{sec2});
    
  X_M_SOE@{sec}@{sec2} = Xss1_M(1,@{sec},@{sec2});
  X_M_ROW@{sec}@{sec2} = Xss1_M(2,@{sec},@{sec2});
  
@#endfor


tau_SOE@{sec}ROW = Tau1(1,@{sec},2);
tau_ROW@{sec}SOE = Tau1(2,@{sec},1);

eps_SOE@{sec}ROW = tau_SOE@{sec}ROW; 
eps_ROW@{sec}SOE = tau_ROW@{sec}SOE;

@#endfor

B_SOE_ss = 0;
B_ROW_ss = 0;
eps_B_SOE_ss = 0;
eps_B_ROW_ss = 0;

end;
resid;
steady; 
check;

tariff_SOE1ROW = [1.21, 1.15];
tariff_SOE3ROW = [1.21, 1.15];
tariff_SOE4ROW = [1.21, 1.15];
tariff_SOE5ROW = [1.21, 1.15];
tariff_SOE6ROW = [1.21, 1.15];

shocks;
var eps_SOE1ROW; 
periods 1:2;
values (tariff_SOE1ROW);

var eps_SOE3ROW; 
periods 1:2;
values (tariff_SOE3ROW);

var eps_SOE4ROW; 
periods 1:2;
values (tariff_SOE4ROW);

var eps_SOE5ROW; 
periods 1:2;
values (tariff_SOE4ROW);

var eps_SOE6ROW; 
periods 1:2;
values (tariff_SOE4ROW);


end;


perfect_foresight_setup(periods=500);
perfect_foresight_solver;
%perfect_foresight_solver(stack_solve_algo=7, solve_algo=9);




