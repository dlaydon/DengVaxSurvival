#pragma once

#include "HeaderAndSwitches.h"

DType Choose_K					(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_Sum_Rho_K			(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_Survive_K			(int TrialArm, int BaselineSeroStatus, int country, int Age, int PhaseSeverity, const Params_Struct &PARAMS, char ModelVariant);
DType Choose_BaseHazMult		(DType K_singlesero,	int TrialArm,	int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_BaseHazMult		(						int TrialArm,	int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_SumBaseHazMult		(int TrialArm,							int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_VacHazMult			(DType K_singlesero,					int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_VacHazMult			(										int BaselineSeroStatus, int country, int Age, int PhaseSeverity, int serotype,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_SumVacHazMult		(DType K_singlesero,					int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Choose_SumVacHazMult		(										int BaselineSeroStatus, int country, int Age, int PhaseSeverity,				const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

void Calc_KplusValues_country_PhaseSeverity	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity);
void Calc_KplusValues_PhaseSeverity			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity); 
void Calc_KplusValues_country					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country);
void Calc_KplusValues_All						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_KplusPrimes_country_PhaseSeverity(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity);
void Calc_KplusPrimes_PhaseSeverity		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity);
void Calc_KplusPrimes_country				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country);
void Calc_KplusPrimes_All					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

void Calc_K_Primes_BS_PhaseSeverity_country_serotype(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus, int PhaseSeverity, int country, int serotype);
void Calc_K_Primes_ALL(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_K_Primes_BS(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus);

void Calc_SumRhoK0_SNeg_Primes_country_PhaseSeverity	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int PhaseSeverity);
void Calc_SumRhoK0_SNeg_Primes_ALL					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

void Calc_ParamSeroPrevs_country				(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

void Calc_RhosFrom_qParams						(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhos_country						(int &country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhos_ALL							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffs_country_BaselineSeroStatus(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int BaselineSeroStatus);
void Calc_SumRhoEffs_country					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country); 
void Calc_SumRhoEffs_BaselineSeroStatus			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus); 
void Calc_SumRhoEffs_ALL						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
void Calc_SumRhoEffs_ALL						(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);

void Calc_SumRhoK_country_PhaseSeverity_PrevInf(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity, int PrevInf); 
void Calc_SumRhoKs_country_PhaseSeverity		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity); 
void Calc_SumRhoKs_PhaseSeverity_PrevInf		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity, int PrevInf);
void Calc_SumRhoKs_country_PrevInf				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PrevInf);
void Calc_SumRhoKs_country						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country);
void Calc_SumRhoKs_PhaseSeverity				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity);
void Calc_SumRhoKs_PrevInf						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PrevInf);
void Calc_SumRhoKs_ALL							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoKs_ALL							(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);

void Calc_SumRhoEffKs_c_BS_PS	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity, int BaselineSeroStatus);
void Calc_SumRhoEffKs_country_PhaseSeverity					(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country, int PhaseSeverity);
void Calc_SumRhoEffKs_PhaseSeverity_BaselineSeroStatus			(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity, int BaselineSeroStatus);
void Calc_SumRhoEffKs_country_BaselineSeroStatus				(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int country, int BaselineSeroStatus);
void Calc_SumRhoEffKs_country									(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int &country);
void Calc_SumRhoEffKs_PhaseSeverity							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int PhaseSeverity);
void Calc_SumRhoEffKs_BaselineSeroStatus						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus);
void Calc_SumRhoEffKs_ALL										(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffKs_ALL										(Params_Struct &CurrentPARAMS, Params_Struct &ProposedPARAMS, const Housekeeping_Struct &HOUSE);

void Calc_BaseHazValues						(int &country, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void Calc_IntBaseHazValues					(int &country, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE);
void Calc_IVH_Values_c_BS_PS		(int country, int BaselineSeroStatus, int PhaseSeverity, Params_Struct &PARAMS, const DATA_struct &DATA, const Housekeeping_Struct &HOUSE); //// Tempting to put this in Initialize_Params, but put it in main as otherwise you interupt following flow: Initialize params, use historical hazards to populate strata sets, calculate IntVacHaz array. . 

bool CanEfficacyBeNegative(int BaselineSeroStatus, const Housekeeping_Struct &HOUSE);
bool IsEfficacyNegative					(int BaselineSeroStatus, int PhaseSeverity, int AgeInYears, int serotype, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE) ;//// really asking IS the initial value of "transient immunity" negative, whether there are age, phase/severity, or serotype, serostatus effects. Not asking whether it can in principle be negative. That is covered in function "CanEfficacyBeNegative" above. 

void Calc_SumRhoEffNegs_c_BS_PS	(int country, int BaselineSeroStatus, int PhaseSeverity, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffNegs_country					(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffNegs_PhaseSeverity			(int PhaseSeverity, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffNegs_BaselineSeroStatus		(int BaselineSeroStatus, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_SumRhoEffNegs_All						(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Calc_AggTransImmunity						(int BaselineSeroStatus, int serotype,	int AgeInYears, int PhaseSeverity,									const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE); 
DType Calc_AggTransImmunity						(int BaselineSeroStatus, int serotype,	int AgeInYears, int PhaseSeverity, const EffNegWane_Option ENW_Option,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Calc_SumRhoAggregateTransImmunity			(int BaselineSeroStatus, int country,	int AgeInYears, int PhaseSeverity, const EffNegWane_Option ENW_Option,	const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType Calc_SumRho_Ks_AggregateTransImmunity		(int BaselineSeroStatus, int country, int AgeInYears, int PhaseSeverity, EffNegWane_Option ENW_Option, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType DetermineAltWaneEffBaseHazMult			(int BaselineSeroStatus, int PhaseSeverity, int AgeInYears, int serotype, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_Hospital_Ks_country					(int country, Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_Hospital_Ks							(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Update_Severe_Ks								(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
void Calc_WaningValues		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus, const int &NoDaysOfFollowUp);
void Calc_AgeVaccEffMults	(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE, int BaselineSeroStatus);
void Calc_Age_HazMults		(Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
