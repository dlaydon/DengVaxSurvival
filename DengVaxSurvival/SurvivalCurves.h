#pragma once

#include "StructureDefs.h"
#include "HeaderAndSwitches.h"
#include "ConsolePrinting.h"


void CalculateAttackRates				(DType **SurvivalTable, DType **AttackRates, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int &iter, const Chains_Struct &CHAINS, int FinalDay, DType Duration);
void CalculateAttackRates_PassiveOnly	(DType **SurvivalTable, DType **AttackRates, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int &iter, const Chains_Struct &CHAINS, int FirstDay, int FinalDay, DType Duration);
void GenerateSurvivalCurves				(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, Survival_Struct &SURVIVE, int iter, const Chains_Struct &CHAINS, Augment_Struct &AUG, bool AllPatients);
void GenerateSurvivalCurves(const DATA_struct &DATA, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE,
	Survival_Struct &SURVIVE, int Iteration_OR_AR_Post_Stat_Index, const Chains_Struct &CHAINS, Augment_Struct &AUG, bool AllPatients,
	//// these arguments were put in for Mean and Modal Posterior SCurves and attack rates. In your main code, you don't set them at all - that is, you call the overloaded function which gives the defaults. 
	bool OutputIndividualSurvivalTables, bool OutputIndividual_Passive_SurvivalTables, DType **** AttackRateContainer,
	string ParamSet_SCurve_RootFilename, string ParamSet_PassiveSCurve_RootFilename, bool Use_Default_AR_Post_Sample_No, bool CleanFirstDay_IndividualCurves);

void CalculateSurvivalCurveOutput		(Survival_Struct &SURVIVE);

