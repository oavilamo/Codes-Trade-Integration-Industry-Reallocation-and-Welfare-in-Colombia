Manual for replicating results "Trade Integration, Industry Reallocation, and Welfare in Colombia"

Files Description: 
Model_***.m: Contain the models to run

1. Model_Benchmark.m: 
	Runs the benchmark model from the paper and stores the initial steady state at Ini_ss.mat. This steady state is used as a starting point in some of the extensions.
	It first runs the step by step calbration, starting from a fully symmetric model (countries and sectors) and then additing heterogeneity at the sector and country level.
	It finally runs the dynare code for the transition path and store the results in oo_Het.mat and M_Het.mat.
	The file stores the parameters, initial and final conditions for running the Dynare code at:
		Parameters: Params_AG.mat, Params_CO.mat, Params_S.mat, Params_I.mat, Params_Omega.mat, Params_lambda.mat
		Initial Steady State: Xss0_CO.mat, Xss0_S.mat, Xss0_K.mat, Tau0.mat, Ice0_0.mat, Ice1_0.mat, Xss0_M.mat
		Final Steady State: Xss1_CO.mat, Xss1_S.mat, Xss1_K.mat, Tau1.mat, Ice0_1.mat, Ice1_1.mat, Xss1_M.mat
	The Dynare File for the Benchmark Model is Model_NS_Bond.mod.

2. Steady_***.m:
	Functions that calculate the Steady State for different specifications of the model (symmetric, asymmetric). 
	The extensions Calib find the stedy state including the specific targets from the data.

3. AlternativePolicies.m
	Runs the code for alternative paths for tariffs and sensitivity. By chaning the variable dum_Acc it is possible to run each of the scenarios.
	dum_Acc = 1: Accelerated. Requires the dynare file Model_NS_BondHom.mod. Saves the results in oo_HetAcc.mat and M_HetAcc.mat 
	dum_Acc = 2: Delayed. Requires the dynare file Model_NS_BondDelay.mod. Saves the results in oo_HetDel.mat and M_HetDel.mat 
	dum_Acc = 3: Transitory. Requires the dynare file Model_NS_BondBack.mod. Saves the results in oo_HetBack.mat and M_HetBack.mat 
	dum_Acc = 4: No Habit. Requires the dynare file Model_NS_NoHab.mod. Saves the results in oo_NoHab.mat and M_NoHab.mat 
	dum_Acc = 5: Low Bond Adjustment Costs. Requires the dynare file Model_NS_NoBAC.mod. Saves the results in oo_NoBAC.mat and M_NoBAC.mat 
	dum_Acc = 6: No Bonds. Requires the dynare file Model_NS_NoBond.mod. Saves the results in oo_NoBond.mat and M_NoBond.mat 

4. Sensitivity (Al):
	Model_NoIO.m: Runs the model with roundabout production (No IO). Requires the dynare file Model_NS_Bond.mod. Saves the results in oo_HetNoIO.mat and M_HetNoIO.mat
	Model_CommonIC.m: Runs the model with common technology for consumption and investment. Requires the dynare file Model_NS_Bond.mod. Saves the results in oo_HetIC.mat and M_HetIC.mat
	Model_Static.m: Runs the model with static export decision. Requires the dynare file Model_NSStatic.mod. Two versions: with and without bonds (FA). Saves the results in oo_HetStatic.mat and M_HetStatic.mat
	Model_TwoSector.m: Runs the model with two sectors. Requires the dynare file Model_NS_Bond2.mod. Saves the results in oo_TwoS.mat and M_TwoS.mat 
	Model_Sectortariffs.m: Runs the sector by sector trade liberalization. Change cathegorical variable dum_sector for each sector (1 (agriculture), 3 (consumption), 5 (industrial), and 6 (capital equipment)). 
		Also modify the dynare file "Model_NS_BondSector.mod" by commenting the lines 1142 to 1153 depending on the sector that is being liberalized. Only keep active the three lines that correspond to the particular sector.
		The code saves the results in oo_sectorAG.mat, oo_sectorMM1.mat, oo_sectorMM3.mat, oo_sectorMM4.mat, M_sectorAG.mat, M_sectorMM1.mat, M_sectorMM3.mat, M_sectorMM4.mat
	Model_HomTariffIni.m: Runs the model that makes initial tariffs equal. Requires the dynare file Model_NSBondHom2.mod. Saves the results in oo_HetHomIni.mat and M_HetHomIni.mat
	Model_HomTariffFin.m: Runs the model that makes final tariffs equal. Requires the dynare file Model_NS_BondHomTariffFin.mod. Saves the results in oo_Hom.mat and M_Hom.mat

5. Dynare files:
	Model_NS_Bond.mod: Benchmark model
	Model_NS_BondHom.mod: Accelerated tariff cut
	Model_NS_BondDelay.mod: Delayed tariff cut
	Model_NS_BondBack.mod: Transitory tariff cut
	Model_NS_NoHab.mod: Model with no Habits
	Model_NS_NoBAC.mod: Model with lower bond adjustment costs
	Model_NS_NoBond.mod: Model with No Bonds
	Model_NS_BondStatic.mod: Model with static export decision


6. Welfare.mat:
	You need to run first all the Models.
	Loads the results from the different scenarios and calculates the static and dynamic welfare changes for each of them. 
	Replicates Table 4.

7. The additional folders are results of running Dynare

*****************************************************************************************
To run the model:

First run the Model_Benchmark.m
Then run AlternativePolicies.m and Sensitivity (Model by Model)

